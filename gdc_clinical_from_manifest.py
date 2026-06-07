#!/usr/bin/env python3
"""
Fetch clinical metadata for exactly the RNA-seq files in a GDC manifest.

Inputs:
  - A GDC manifest file (TSV) that contains a column 'id' (UUIDs). 
    Works even if the column is called 'file_id'.

Outputs (in outdir):
  - file_to_case_map.csv       : map of file_id -> case_id, submitter_id, project_id, sample submitter_id
  - clinical_cases_flat.csv    : one row per case with commonly used clinical fields
  - clinical_cases_raw.jsonl   : raw per-case JSON (one line per case)
  - missing_file_ids.txt       : file UUIDs that couldn’t be resolved (if any)

Usage:
  python gdc_clinical_from_manifest.py --manifest manifest.txt --outdir PAAD_Clinical
"""

import argparse
import json
import math
import os
import sys
from typing import List, Dict, Any, Tuple, Iterable

import pandas as pd
import requests
from tqdm import tqdm

GDC_BASE = "https://api.gdc.cancer.gov"
SESSION = requests.Session()
SESSION.headers.update({"Content-Type": "application/json"})

def chunked(iterable: Iterable, n: int) -> Iterable[List]:
    chunk = []
    for x in iterable:
        chunk.append(x)
        if len(chunk) == n:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def read_manifest_file_ids(path: str) -> List[str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    col_candidates = [c for c in df.columns if c.lower() in {"id", "file_id"}]
    if not col_candidates:
        raise ValueError(
            f"Manifest at {path} lacks an 'id' or 'file_id' column. "
            "Download a proper manifest from the GDC Portal."
        )
    file_id_col = col_candidates[0]
    ids = df[file_id_col].dropna().unique().tolist()
    if not ids:
        raise ValueError("No file IDs found in manifest.")
    return ids

def gdc_post(endpoint: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    url = f"{GDC_BASE}/{endpoint.lstrip('/')}"
    r = SESSION.post(url, data=json.dumps(payload), timeout=120)
    r.raise_for_status()
    return r.json()

def fetch_files_mapping(file_ids: List[str]) -> pd.DataFrame:
    """
    Return rows with: file_id, case_id, case_submitter_id, project_id, sample_submitter_id
    """
    fields = [
        "file_id",
        "cases.case_id",
        "cases.submitter_id",
        "cases.project.project_id",
        "cases.samples.sample_id",
        "cases.samples.sample_type",
        "cases.samples.submitter_id",
    ]
    rows = []
    missing = []

    # GDC Filters support up to a lot, but be gentle (e.g., 200 IDs per batch)
    for batch in tqdm(list(chunked(file_ids, 200)), desc="Resolving file->case"):
        payload = {
            "filters": {
                "op": "in",
                "content": {"field": "file_id", "value": batch},
            },
            "fields": ",".join(fields),
            "format": "JSON",
            "size": 5000,
        }
        data = gdc_post("/files", payload)
        hits = data.get("data", {}).get("hits", [])
        seen = set()
        for h in hits:
            fid = h.get("file_id")
            seen.add(fid)
            cases = h.get("cases", []) or []
            if not cases:
                rows.append({
                    "file_id": fid,
                    "case_id": None,
                    "case_submitter_id": None,
                    "project_id": None,
                    "sample_submitter_id": None,
                    "sample_type": None,
                })
                continue
            # Usually one case per file; still iterate defensively
            for c in cases:
                samples = c.get("samples", []) or [dict()]
                if not samples:
                    samples = [dict()]
                for s in samples:
                    rows.append({
                        "file_id": fid,
                        "case_id": c.get("case_id"),
                        "case_submitter_id": c.get("submitter_id"),
                        "project_id": (c.get("project") or {}).get("project_id"),
                        "sample_submitter_id": s.get("submitter_id"),
                        "sample_type": s.get("sample_type"),
                    })
        # flag any not returned
        for fid in batch:
            if fid not in seen:
                missing.append(fid)

    df = pd.DataFrame(rows).drop_duplicates()
    return df, missing

def fetch_clinical_for_cases(case_ids: List[str]) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
    """
    Pull rich clinical fields via /cases with expansions and flatten a tidy table.
    """
    # request common clinical fields; add more if needed
    fields = [
        "case_id",
        "submitter_id",
        "project.project_id",
        # demographics
        "demographic.ethnicity",
        "demographic.gender",
        "demographic.race",
        "demographic.year_of_birth",
        "demographic.year_of_death",
        # diagnoses (first/primary will be flattened preferentially)
        "diagnoses.age_at_diagnosis",
        "diagnoses.tumor_stage",
        "diagnoses.classification_of_tumor",
        "diagnoses.primary_diagnosis",
        "diagnoses.tumor_grade",
        "diagnoses.morphology",
        "diagnoses.site_of_resection_or_biopsy",
        "diagnoses.days_to_death",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.vital_status",
        "diagnoses.prior_malignancy",
        "diagnoses.prior_treatment",
        "diagnoses.residual_disease",
        "diagnoses.treatments.therapy_type",
        "diagnoses.treatments.treatment_type",
        "diagnoses.treatments.treatment_anatomic_site",
        "diagnoses.treatments.treatment_effect",
        "diagnoses.treatments.treatment_outcome",
        # exposures (sometimes empty)
        "exposures.alcohol_history",
        "exposures.bmi",
        "exposures.height",
        "exposures.weight",
        "exposures.tobacco_smoking_history",
    ]
    raw_cases = []
    flat_rows = []

    # batch case_ids to avoid very long filters
    for batch in tqdm(list(chunked(case_ids, 200)), desc="Fetching clinical (/cases)"):
        payload = {
            "filters": {
                "op": "in",
                "content": {"field": "case_id", "value": batch},
            },
            "fields": ",".join(fields),
            "expand": "diagnoses,diagnoses.treatments,demographic,exposures",
            "format": "JSON",
            "size": 5000,
        }
        data = gdc_post("/cases", payload)
        hits = data.get("data", {}).get("hits", [])
        for h in hits:
            raw_cases.append(h)  # keep intact for JSONL

            # flatten: pick a "primary" diagnosis (first if multiple)
            diags = h.get("diagnoses") or []
            d0 = diags[0] if diags else {}

            # Optionally aggregate treatments into a semicolon list
            trts = d0.get("treatments") or []
            def join_nonnull(seq, key):
                vals = [str(x.get(key)) for x in seq if x and x.get(key) not in [None, "", "None"]]
                return ";".join(vals) if vals else None

            flat_rows.append({
                "case_id": h.get("case_id"),
                "case_submitter_id": h.get("submitter_id"),
                "project_id": (h.get("project") or {}).get("project_id"),

                # demographics
                "gender": (h.get("demographic") or {}).get("gender"),
                "race": (h.get("demographic") or {}).get("race"),
                "ethnicity": (h.get("demographic") or {}).get("ethnicity"),
                "year_of_birth": (h.get("demographic") or {}).get("year_of_birth"),
                "year_of_death": (h.get("demographic") or {}).get("year_of_death"),

                # diagnoses (primary)
                "primary_diagnosis": d0.get("primary_diagnosis"),
                "tumor_stage": d0.get("tumor_stage"),
                "tumor_grade": d0.get("tumor_grade"),
                "classification_of_tumor": d0.get("classification_of_tumor"),
                "morphology": d0.get("morphology"),
                "site_of_resection_or_biopsy": d0.get("site_of_resection_or_biopsy"),
                "age_at_diagnosis_days": d0.get("age_at_diagnosis"),
                "vital_status": d0.get("vital_status"),
                "days_to_death": d0.get("days_to_death"),
                "days_to_last_follow_up": d0.get("days_to_last_follow_up"),
                "prior_malignancy": d0.get("prior_malignancy"),
                "prior_treatment": d0.get("prior_treatment"),
                "residual_disease": d0.get("residual_disease"),

                # treatments (aggregated)
                "treatment_type_list": join_nonnull(trts, "treatment_type"),
                "therapy_type_list": join_nonnull(trts, "therapy_type"),
                "treatment_anatomic_site_list": join_nonnull(trts, "treatment_anatomic_site"),
                "treatment_effect_list": join_nonnull(trts, "treatment_effect"),
                "treatment_outcome_list": join_nonnull(trts, "treatment_outcome"),

                # exposures (first if multiple)
                "bmi": (h.get("exposures") or [{}])[0].get("bmi") if h.get("exposures") else None,
                "alcohol_history": (h.get("exposures") or [{}])[0].get("alcohol_history") if h.get("exposures") else None,
                "height": (h.get("exposures") or [{}])[0].get("height") if h.get("exposures") else None,
                "weight": (h.get("exposures") or [{}])[0].get("weight") if h.get("exposures") else None,
                "tobacco_smoking_history": (h.get("exposures") or [{}])[0].get("tobacco_smoking_history") if h.get("exposures") else None,
            })

    flat_df = pd.DataFrame(flat_rows).drop_duplicates()
    return raw_cases, flat_df

def main():
    ap = argparse.ArgumentParser(description="Fetch GDC clinical metadata for files in a manifest.")
    ap.add_argument("--manifest", required=True, help="Path to GDC manifest (TSV) with 'id' column (file UUIDs).")
    ap.add_argument("--outdir", required=True, help="Output directory, e.g., PAAD_Clinical")
    ap.add_argument("--timeout", type=int, default=120, help="HTTP timeout seconds")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # read file UUIDs
    file_ids = read_manifest_file_ids(args.manifest)
    print(f"Loaded {len(file_ids)} file IDs from manifest")

    # Map file -> case
    file_case_df, missing = fetch_files_mapping(file_ids)
    file_case_path = os.path.join(args.outdir, "file_to_case_map.csv")
    file_case_df.to_csv(file_case_path, index=False)
    print(f"Saved file-to-case map: {file_case_path}")

    if missing:
        miss_path = os.path.join(args.outdir, "missing_file_ids.txt")
        with open(miss_path, "w") as f:
            f.write("\n".join(missing))
        print(f"WARNING: {len(missing)} file IDs not resolved. See: {miss_path}")

    # Unique case_ids for which we want clinical metadata
    case_ids = sorted(x for x in file_case_df["case_id"].dropna().unique().tolist())
    if not case_ids:
        print("No case_ids resolved; nothing to fetch. Check your manifest or network.")
        sys.exit(1)
    print(f"Found {len(case_ids)} unique cases")

    # Fetch clinical metadata and flatten
    raw_cases, flat_df = fetch_clinical_for_cases(case_ids)

    # Save raw JSONL
    raw_path = os.path.join(args.outdir, "clinical_cases_raw.jsonl")
    with open(raw_path, "w") as f:
        for rc in raw_cases:
            f.write(json.dumps(rc, separators=(",", ":"), ensure_ascii=False) + "\n")
    print(f"Saved raw clinical JSONL: {raw_path}")

    # Save flat CSV
    flat_path = os.path.join(args.outdir, "clinical_cases_flat.csv")
    flat_df.to_csv(flat_path, index=False)
    print(f"Saved flat clinical CSV: {flat_path}")

    print("\nDone.")

if __name__ == "__main__":
    try:
        main()
    except requests.HTTPError as e:
        print(f"HTTP error from GDC API: {e} | response: {getattr(e, 'response', None) and e.response.text}")
        sys.exit(2)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(3)
