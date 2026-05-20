# ==============================
# FULL IMMUNE NETWORK PIPELINE
# ==============================

# -------------------------------
# 0) Load libraries
# -------------------------------
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor
BiocManager::install(c("clusterProfiler","org.Hs.eg.db","AnnotationDbi","RCy3","EnsDb.Hsapiens.v86"), ask = FALSE)

# CRAN
install.packages(c("data.table","dplyr","igraph","httr","jsonlite","foreach","doParallel"))

# Load libraries
library(data.table)
library(dplyr)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RCy3)
library(EnsDb.Hsapiens.v86)
library(httr)
library(jsonlite)
library(foreach)
library(doParallel)

# -------------------------------
# 1) Input / Output paths
# -------------------------------
# DESeq2 table
deg_deseq_csv <- "A:/thesis/network_building_2/Combined_DEGs_all_sorted_deseq2.csv"  

# Limma table
deg_limma_csv <- "A:/thesis/network_building_2/MASTER_all_significant_DEGs_limma.csv" 

# Expression matrix (rows = genes, cols = samples)
expr_csv <- "A:/thesis/network_building_2/DESeq2_normalized_counts_all_genes.csv"  

# Immune gene list (optional)
immune_gene_csv <- "A:/thesis/network_building_2/Merged_Immune_Genes.txt"  

# CIBERSORTx output
cibersort_table <- "A:/thesis/network_building_2/CIBERSORTx_Job1_Results.csv"  

# Output folder
output_dir <- "A:/thesis/network_building_2/network_pipeline_outputs"
if(!dir.exists(output_dir)) dir.create(output_dir)

# -------------------------------
# 2) Thresholds
# -------------------------------
cor_threshold <- 0.6
pval_threshold <- 0.05
string_score_cut <- 700
species_taxon <- 9606 # human

# -------------------------------
# 3) Load DE tables & filter DEGs
# -------------------------------

# DESeq2
deg_raw <- fread(deg_deseq_csv)

deg_deseq <- deg_raw[padj < 0.05 & abs(log2FoldChange) >= 0.5]
cat("DESeq2 filtered DEGs:", nrow(deg_deseq), "\n")

# Limma
deg_raw2 <- fread(deg_limma_csv)
deg_limma <- deg_raw2[adj.P.Val < 0.05 & abs(logFC) >= 0.5]
cat("Limma filtered DEGs:", nrow(deg_limma), "\n")

# -------------------------------
# 4) Candidate genes: union + intersection
# -------------------------------
sig_deseq <- toupper(deg_deseq$gene)
sig_limma <- toupper(deg_limma$gene)

set_intersection <- intersect(sig_deseq, sig_limma)
set_union <- union(sig_deseq, sig_limma)

cat("Intersection DEGs:", length(set_intersection), "\n")
cat("Union DEGs:", length(set_union), "\n")

# Optional: load immune genes
immune_genes <- if(file.exists(immune_gene_csv)){
  unique(trimws(scan(immune_gene_csv, what=character(), sep="\n", quiet=TRUE)))
} else character(0)

candidates <- unique(c(set_union, toupper(immune_genes)))
candidates <- candidates[candidates != ""]
cat("Total candidate genes:", length(candidates), "\n")

library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(dplyr)

# candidates = vector of gene symbols (uppercase)
candidates <- toupper(candidates)

# ---- 1) org.Hs.eg.db mapping ----
# Only SYMBOL, ENTREZID, GENENAME are valid
ann1 <- AnnotationDbi::select(
  x = org.Hs.eg.db,
  keys = candidates,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "SYMBOL"
)
ann1 <- ann1[!duplicated(ann1$SYMBOL), ]

# ---- 2) EnsDb mapping ----
# Use SYMBOL -> GENEID + GENEBIOTYPE
ann_ens <- AnnotationDbi::select(
  x = EnsDb.Hsapiens.v86,
  keys = candidates,
  columns = c("SYMBOL", "GENEID", "GENEBIOTYPE"),
  keytype = "SYMBOL"
)
ann_ens <- ann_ens[!duplicated(ann_ens$SYMBOL), ]

# ---- 3) Merge safely ----
ann <- left_join(ann1, ann_ens, by = "SYMBOL")

# rename for clarity
colnames(ann)[colnames(ann) == "GENEID"] <- "ensembl_gene_id"
colnames(ann)[colnames(ann) == "GENEBIOTYPE"] <- "biotype"

head(ann)



protein_coding_symbols <- unique(ann$SYMBOL[ann$biotype == "protein_coding"])
cat("Protein-coding genes:", length(protein_coding_symbols), "\n")

# Intersect with DEGs
sig_deseq_pc <- intersect(sig_deseq, protein_coding_symbols)
sig_limma_pc <- intersect(sig_limma, protein_coding_symbols)
inter_pc <- intersect(sig_deseq_pc, sig_limma_pc)
union_pc <- union(sig_deseq_pc, sig_limma_pc)

# Corrected DEGs (choose conservative DESeq2)
corrected_degs <- if(length(sig_deseq_pc) >= 50) sig_deseq_pc else if(length(inter_pc) >= 50) inter_pc else union_pc
cat("Corrected DEGs used:", length(corrected_degs), "\n")
fwrite(data.table(gene=corrected_degs), file.path(output_dir, "corrected_DEGs.txt"))

# -------------------------------
# 6) Load expression matrix
# -------------------------------
expr_dt <- fread(expr_csv)
expr_mat <- as.matrix(expr_dt[, -1, with=FALSE])
rownames(expr_mat) <- toupper(expr_dt[[1]])

# -------------------------------
# 7) Load CIBERSORTx
# -------------------------------
#cib_dt <- fread(cibersort_table)
#idcol <- intersect(names(cib_dt), c("Sample","sample","Mixture","mixture","ID","id"))
#if(length(idcol) >= 1) {
#  sidcol <- idcol[1]
 # rownames(cib_dt) <- cib_dt[[sidcol]]
  #cib_mat <- as.matrix(cib_dt[, setdiff(names(cib_dt), sidcol), with=FALSE])
#} else {
#  rownames(cib_dt) <- cib_dt[[1]]
#  cib_mat <- as.matrix(cib_dt[, -1, with=FALSE])
#}

# Load CIBERSORTx output
cib_dt <- fread(cibersort_table)
# First column is sample IDs
sample_col <- names(cib_dt)[1]
sample_ids <- cib_dt[[sample_col]]  # save sample IDs

# Keep only immune cell columns (numeric)
immune_cols <- setdiff(names(cib_dt), sample_col)
cib_mat <- as.matrix(cib_dt[, ..immune_cols, with=FALSE])

# Assign rownames as sample IDs
rownames(cib_mat) <- sample_ids

# Align samples
common_samps <- intersect(colnames(expr_mat), rownames(cib_mat))
expr_mat <- expr_mat[, common_samps]
cib_mat <- cib_mat[common_samps, ]

# -------------------------------
# # -------------------------------
# 8) Master gene list (DEGs + immune)
# -------------------------------
master_genes_fixed <- intersect(
  unique(c(corrected_degs, toupper(immune_genes))), 
  rownames(expr_mat)
)
cat("Master genes in expression matrix:", length(master_genes_fixed), "\n")
fwrite(data.table(gene=master_genes_fixed), file.path(output_dir, "master_genes_fixed.txt"))

# Limit number of master genes
max_keep <- 600
if(length(master_genes_fixed) > max_keep){
  gene_vars <- apply(expr_mat[master_genes_fixed, , drop=FALSE], 1, var)
  ord <- order(gene_vars, decreasing=TRUE)
  master_genes <- master_genes_fixed[ord][1:max_keep]
} else {
  master_genes <- master_genes_fixed
}
cat("Final master genes used for correlation:", length(master_genes), "\n")

# -------------------------------
# 9) Co-expression (Pearson correlation)
# -------------------------------
expr_mat_num <- as.matrix(expr_mat)
mode(expr_mat_num) <- "numeric"
rownames(expr_mat_num) <- rownames(expr_mat)

expr_sub <- expr_mat_num[master_genes, , drop=FALSE]
expr_t <- t(expr_sub)
n_samples <- nrow(expr_t)

cor_mat <- cor(expr_t, method="pearson", use="pairwise.complete.obs")

library(data.table)
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("geneA", "geneB", "cor"))
cor_dt <- cor_dt[geneA < geneB]  # upper triangle only

df <- n_samples - 2
cor_dt[, tval := cor * sqrt(df / (1 - cor^2))]
cor_dt[, pval := 2 * pt(-abs(tval), df=df)]
cor_dt[, tval := NULL]

fwrite(cor_dt, file.path(output_dir, "correlation_long.csv"))
cat("Correlation calculation complete. Total pairs:", nrow(cor_dt), "\n")

# -------------------------------
# 10) Filter correlations
# -------------------------------
cor_dt[, p_adj := p.adjust(pval, method="BH")]
rcut <- 0.7
key_genes <- unique(c(corrected_degs, intersect(toupper(immune_genes), rownames(expr_sub))))
cor_final <- cor_dt[p_adj < 0.05 & abs(cor) >= rcut & (geneA %in% key_genes | geneB %in% key_genes)]

max_edges <- 20000
if(nrow(cor_final) > max_edges){
  cor_final <- cor_final[order(-abs(cor))][1:max_edges]
}

fwrite(cor_final[, .(source=geneA, target=geneB, cor, pval, p_adj)],
       file.path(output_dir, "expr_expr_edges_filtered.csv"))

# -------------------------------
# 11) Gene <-> Immune correlations
# -------------------------------
genes_for_immune_corr <- intersect(key_genes, rownames(expr_sub))
gene_imm_list <- list()
for(g in genes_for_immune_corr){
  x <- as.numeric(expr_sub[g, ])
  for(cell in colnames(cib_mat)){
    y <- as.numeric(cib_mat[, cell])
    r <- cor(x, y)
    tval <- r * sqrt((length(x)-2)/(1-r^2))
    pval <- 2*pt(-abs(tval), df=length(x)-2)
    gene_imm_list[[length(gene_imm_list)+1]] <- data.table(gene=g, immune_cell=cell, cor=r, pval=pval)
  }
}
gene_imm_dt <- rbindlist(gene_imm_list)
gene_imm_dt[, p_adj := p.adjust(pval, method="BH")]
gene_imm_edges <- gene_imm_dt[abs(cor) >= 0.6 & p_adj < 0.05]
fwrite(gene_imm_edges, file.path(output_dir, "gene_immune_edges_filtered.csv"))

# -------------------------------
# 12) STRING interactions
# -------------------------------
fetch_string <- function(genes_vector, species=9606, required_score=700){
  chunks <- split(genes_vector, ceiling(seq_along(genes_vector)/200))
  all_res <- list()
  for(i in seq_along(chunks)){
    ids <- paste(chunks[[i]], collapse="%0D")
    url <- sprintf("https://string-db.org/api/tsv/network?identifiers=%s&species=%s&required_score=%s",
                   ids, species, required_score)
    resp <- httr::GET(url)
    if(resp$status_code != 200) next
    txt <- httr::content(resp, "text", encoding="UTF-8")
    dt <- tryCatch(fread(txt), error=function(e) NULL)
    if(!is.null(dt)) all_res[[length(all_res)+1]] <- dt
    Sys.sleep(1)
  }
  if(length(all_res) == 0) return(data.table())
  return(rbindlist(all_res, fill=TRUE))
}

string_dt <- fetch_string(unique(toupper(master_genes_fixed)), species=species_taxon, required_score=string_score_cut)
if(nrow(string_dt) > 0){
  setnames(string_dt, tolower(names(string_dt)))
  string_edges <- string_dt[, .(protein1=preferredname_a, protein2=preferredname_b, score=score)]
  string_edges[, score01 := score/1000]
  fwrite(string_edges, file.path(output_dir, "string_interactions.tsv"))
}

# -------------------------------
# 13) Build network for Cytoscape
# -------------------------------
edge_list <- list()
if(exists("cor_final")) edge_list[[length(edge_list)+1]] <- cor_final[, .(source=geneA, target=geneB, interaction="expr_corr", score=abs(cor), pval, p_adj)]
if(exists("gene_imm_edges")) edge_list[[length(edge_list)+1]] <- gene_imm_edges[, .(source=gene, target=immune_cell, interaction="immune_corr", score=abs(cor), pval, p_adj)]
if(exists("string_edges")) edge_list[[length(edge_list)+1]] <- string_edges[, .(source=protein1, target=protein2, interaction="STRING", score=score01, pval=NA, p_adj=NA)]
edges_all <- rbindlist(edge_list, fill=TRUE)
fwrite(edges_all, file.path(output_dir, "network_edges_for_cytoscape.csv"))

node_ids <- unique(c(edges_all$source, edges_all$target))
nodes <- data.table(id=node_ids)
nodes[, in_DEG := toupper(id) %in% toupper(corrected_degs)]
nodes[, in_immune_list := toupper(id) %in% toupper(immune_genes)]
nodes[, type := fifelse(in_immune_list, "immune_gene", fifelse(in_DEG, "DEG", "other"))]
nodes[, avg_expr := {
  idx <- match(toupper(id), toupper(rownames(expr_mat)))
  ifelse(!is.na(idx), rowMeans(expr_mat[idx,,drop=FALSE]), NA_real_)
}]
fwrite(nodes, file.path(output_dir, "network_nodes_for_cytoscape.csv"))

# -------------------------------
# 14) Centrality and hub genes
# -------------------------------
library(igraph)
g <- graph_from_data_frame(d=edges_all[, .(source, target)], directed=FALSE, vertices=nodes)
deg <- degree(g)
bet <- betweenness(g)
clo <- closeness(g, mode="all", normalized=TRUE)
node_metrics <- data.table(id=names(deg), degree=deg, betweenness=bet, closeness=clo)
nodes_metrics <- merge(nodes, node_metrics, by="id", all.x=TRUE)
fwrite(nodes_metrics, file.path(output_dir, "node_metrics.csv"))

deg_cut <- quantile(node_metrics$degree, 0.95, na.rm=TRUE)
bet_cut <- quantile(node_metrics$betweenness, 0.90, na.rm=TRUE)
hubs <- node_metrics[degree >= deg_cut & betweenness >= bet_cut]
fwrite(hubs, file.path(output_dir, "candidate_hub_genes.csv"))
cat("Hub genes identified:", nrow(hubs), "\n")

# -------------------------------
# 15) Save session info
# -------------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Pipeline finished. All outputs are in:", output_dir, "\n")
