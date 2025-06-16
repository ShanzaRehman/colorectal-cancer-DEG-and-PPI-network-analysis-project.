#Step 1: Data Download and Loading

# Load required libraries
library(GEOquery)
library(DESeq2)
library(limma)
library(edgeR)
library(tidyverse)
library(biomaRt)
library(dplyr)
library(stringr)
library(STRINGdb)

library(AnnotationDbi)
library(org.Hs.eg.db)

gse <- getGEO("GSE50760", GSEMatrix = TRUE, getGPL = FALSE)
pheno_data <- pData(gse[[1]])
expr_data <- read_tsv("D:/R programs/GSE50760_raw_counts_GRCh38.p13_NCBI.tsv")

cat("Expression matrix dimensions:", dim(expr_data), "\n")
cat("Phenotype data dimensions:", dim(pheno_data), "\n")


# Find samples in expr_data not in pheno_data
mismatch <- setdiff(colnames(expr_data), rownames(pheno_data))
print(mismatch)

# Remove the unmatched column from expr_data
expr_data_clean <- expr_data[, colnames(expr_data) %in% rownames(pheno_data)]

# Reorder columns of expr_data_clean to match the rows of pheno_data
expr_data_clean <- expr_data_clean[, rownames(pheno_data)]

cat("Cleaned expression matrix dimensions:", dim(expr_data_clean), "\n")
cat("Phenotype data dimensions:", dim(pheno_data), "\n")

# Preview data
head(expr_data_clean[, 1:5])
head(pheno_data[, c("title", "characteristics_ch1")])

pheno_data_df <- as.data.frame(pheno_data)


#Step 2: Data Preprocessing and Quality Control
# Extract sample information
sample_info <- pheno_data_df %>%
  dplyr::select(title, characteristics_ch1) %>%
  dplyr::mutate(
    sample_type = ifelse(grepl("normal|Normal", characteristics_ch1), 
                         "Normal", "Tumor"),
    patient_id = str_extract(title, "\\d+")
  )

# Create design matrix
design_matrix <- data.frame(
  sample_id = rownames(sample_info),
  condition = factor(sample_info$sample_type, levels = c("Normal", "Tumor")),
  patient = sample_info$patient_id
)


# Quality control plots
library(ggplot2)
library(pheatmap)

# 1. Distribution of expression values
expr_long <- expr_data %>%
  as.data.frame() %>%
  rownames_to_column("probe") %>%
  pivot_longer(-probe, names_to = "sample", values_to = "expression")


# Expression distribution plot
p1 <- ggplot(expr_long, aes(x = expression, color = sample)) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Expression Distribution Across Samples") +
  theme(legend.position = "none")

# 2. Sample correlation heatmap
cor_matrix <- cor(expr_data_clean, use = "complete.obs")
pheatmap(cor_matrix, 
         annotation_col = data.frame(
           Condition = design_matrix$condition,
           row.names = colnames(expr_data_clean)
         ),
         main = "Sample Correlation Heatmap")


# 3. Principal Component Analysis
pca_data <- prcomp(t(expr_data_clean), scale = TRUE)

# Remove genes with zero variance
non_zero_var_genes <- apply(expr_data_clean, 1, function(x) var(x) != 0)
expr_data_pca_ready <- expr_data_clean[non_zero_var_genes, ]

# Perform PCA
pca_data <- prcomp(t(expr_data_pca_ready), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  Condition = design_matrix$condition
)

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Samples",
       x = paste0("PC1 (", round(summary(pca_data)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_data)$importance[2,2]*100, 1), "%)"))


print(p1)
print(p2)

#Step 3: Differential Expression Analysis

# Convert expression data to integer counts (simulate RNA-seq data)

# In real RNA-seq, you would start with raw count data
#count_data <- round(2^expr_data_pca_ready)

counts <- read.delim("D:/R programs/GSE50760_raw_counts_GRCh38.p13_NCBI.tsv", 
                     header = TRUE, 
                     row.names = 1, 
                     check.names = FALSE, 
                     sep = "\t")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = design_matrix,
  design = ~ condition
)

#Filter low count genes
keep<- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#RUN DESeq2 analysis
dds<- DESeq(dds)

#Extract resuts
deseq2_results <- results(dds, contrast = c("condition", "Tumor", "Normal"))
deseq2_results <- deseq2_results[complete.cases(deseq2_results),]

#Add gene annotations
deseq2_df <- as.data.frame(deseq2_results) %>%
  rownames_to_column("probe_id") %>%
  mutate(
    gene_symbol= mapIds(org.Hs.eg.db, keys = probe_id,
                        column = "SYMBOL", keytype = "ENTREZID",
                        multiVals = "first")
  )%>%
  arrange(padj)

#Significant genes (DESeq2)
deseq2_sig <- deseq2_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1)

cat("DESeq2: Sifnificant DEGS:", nrow(deseq2_sig), "\n")
    
  
#3.2 limma Analysis
# Design matrix for limma
design_limma <- model.matrix(~ condition, data = design_matrix)

# Fit linear model
fit <- lmFit(expr_data_pca_ready, design_limma)
fit <- eBayes(fit)

# Extract results
limma_results <- topTable(fit, coef = "conditionTumor", 
                          number = Inf, adjust.method = "BH")

# Add fold change and significance flags
limma_df <- limma_results %>%
  rownames_to_column("probe_id") %>%
  mutate(
    log2FoldChange = logFC,
    padj = adj.P.Val,
    gene_symbol = mapIds(
      org.Hs.eg.db,
      keys = probe_id,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
  ) %>%
  arrange(padj)

# Significant genes (limma)
limma_sig <- limma_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("limma: Significant DEGs:", nrow(limma_sig), "\n")

# 3 edgeR Analysis
# Create DGEList object
count_data <- read.delim("D:/R programs/GSE50760_raw_counts_GRCh38.p13_NCBI.tsv", 
                                   header = TRUE, 
                                   row.names = 1, 
                                   check.names = FALSE, 
                                   sep = "\t")
dge <- DGEList(counts = count_data, group = design_matrix$condition)

# Filter low-expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Estimate dispersions
design_edger <- model.matrix(~ condition, data = design_matrix)
dge <- estimateDisp(dge, design_edger)

# Fit quasi-likelihood model
fit_edger <- glmQLFit(dge, design_edger)
qlf <- glmQLFTest(fit_edger, coef = "conditionTumor")

# Extract results
edger_results <- topTags(qlf, n = Inf)$table

# Format results
edger_df <- edger_results %>%
  rownames_to_column("probe_id") %>%
  mutate(
    padj = FDR,
    log2FoldChange = logFC,
    gene_symbol = mapIds(org.Hs.eg.db, keys = probe_id, 
                         column = "SYMBOL", keytype = "ENTREZID", 
                         multiVals = "first")
  ) %>%
  arrange(padj)

# Significant genes (edgeR)
edger_sig <- edger_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

cat("edgeR: Significant DEGs:", nrow(edger_sig), "\n")

#Step 4: Comparison of Methods and Visualization
# Compare results across methods
library(VennDiagram)
library(grid)

# Get gene lists
deseq2_genes <- deseq2_sig$gene_symbol[!is.na(deseq2_sig$gene_symbol)]
limma_genes <- limma_sig$gene_symbol[!is.na(limma_sig$gene_symbol)]
edger_genes <- edger_sig$gene_symbol[!is.na(edger_sig$gene_symbol)]

# Venn diagram
venn.plot <- venn.diagram(
  list(DESeq2 = deseq2_genes, limma = limma_genes, edgeR = edger_genes),
  filename = NULL,
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontfamily = "serif",
  main = "Overlap of DEGs between Methods"
)
grid.newpage()      # Start a new plot page
grid.draw(venn.plot)  # Draw the Venn diagram

# Find consensus genes (detected by at least 2 methods)
all_genes <- union(union(deseq2_genes, limma_genes), edger_genes)
consensus_genes <- c()

for(gene in all_genes) {
  count <- sum(c(gene %in% deseq2_genes, 
                 gene %in% limma_genes, 
                 gene %in% edger_genes))
  if(count >= 2) {
    consensus_genes <- c(consensus_genes, gene)
  }
}

cat("Consensus DEGs (detected by ≥2 methods):", length(consensus_genes), "\n")


# Volcano plots
create_volcano_plot <- function(data, title, fc_col, pval_col) {
  data$significance <- "Not Significant"
  data$significance[data[[pval_col]] < 0.05 & abs(data[[fc_col]]) > 1] <- "Significant"
  
  ggplot(data, aes_string(x = fc_col, y = paste0("-log10(", pval_col, ")"), 
                          color = "significance")) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)
}

# Create volcano plots
p_volcano_deseq2 <- create_volcano_plot(deseq2_df, "DESeq2 Volcano Plot", 
                                        "log2FoldChange", "padj")
p_volcano_limma <- create_volcano_plot(limma_df, "limma Volcano Plot", 
                                       "log2FoldChange", "padj")
p_volcano_edger <- create_volcano_plot(edger_df, "edgeR Volcano Plot", 
                                       "log2FoldChange", "padj")

print(p_volcano_deseq2)
print(p_volcano_limma)
print(p_volcano_edger)

#Step 5: Functional Annotation and Pathway Analysis
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(enrichplot)

# Convert gene symbols to Entrez IDs for pathway analysis
consensus_entrez <- bitr(consensus_genes, fromType = "SYMBOL", 
                         toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Gene Ontology enrichment analysis
go_bp <- enrichGO(gene = consensus_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

go_mf <- enrichGO(gene = consensus_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

go_cc <- enrichGO(gene = consensus_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)


# KEGG pathway analysis
kegg_pathways <- enrichKEGG(gene = consensus_entrez$ENTREZID,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)

# Reactome pathway analysis
reactome_pathways <- enrichPathway(gene = consensus_entrez$ENTREZID,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   readable = TRUE)

# Visualization of enrichment results
# GO Biological Process
p_go_bp <- dotplot(go_bp, showCategory = 20) + 
  ggtitle("GO Biological Process Enrichment")

# KEGG pathways
p_kegg <- dotplot(kegg_pathways, showCategory = 15) + 
  ggtitle("KEGG Pathway Enrichment")

# Reactome pathways
p_reactome <- dotplot(reactome_pathways, showCategory = 15) + 
  ggtitle("Reactome Pathway Enrichment")

print(p_go_bp)
print(p_kegg)
print(p_reactome)

# Gene-concept network
if(nrow(as.data.frame(kegg_pathways)) > 0) {
  p_cnet <- cnetplot(kegg_pathways, categorySize="pvalue", foldChange = NULL)
  print(p_cnet)
}
#Step 6: Protein-Protein Interaction (PPI) Network Analysis
library(STRINGdb)

# Initialize STRING database
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400)

# Map genes to STRING IDs
consensus_mapped <- string_db$map(data.frame(gene = consensus_genes), 
                                  "gene", removeUnmappedRows = TRUE)

# Get PPI network
ppi_network <- string_db$get_interactions(consensus_mapped$STRING_id)

# Network statistics
cat("PPI Network Statistics:\n")
cat("Nodes:", length(unique(c(ppi_network$from, ppi_network$to))), "\n")
cat("Edges:", nrow(ppi_network), "\n")

# Find hub genes (genes with many interactions)
node_degrees <- table(c(ppi_network$from, ppi_network$to))
hub_genes <- names(sort(node_degrees, decreasing = TRUE))[1:10]

# Map back to gene symbols
hub_gene_symbols <- consensus_mapped[consensus_mapped$STRING_id %in% hub_genes, c("gene", "STRING_id")]
print("Top 10 Hub Genes:")
print(hub_gene_symbols)

# Plot PPI network (simplified)
if(nrow(ppi_network) > 0) {
  string_db$plot_network(consensus_mapped$STRING_id[1:50])  # Plot subset for clarity
}
#Step 7: Results Summary and Export
# Create comprehensive results summary
results_summary <- data.frame(
  Method = c("DESeq2", "limma", "edgeR", "Consensus"),
  Total_DEGs = c(nrow(deseq2_sig), nrow(limma_sig), 
                 nrow(edger_sig), length(consensus_genes)),
  Upregulated = c(
    sum(deseq2_sig$log2FoldChange > 0),
    sum(limma_sig$log2FoldChange > 0),
    sum(edger_sig$log2FoldChange > 0),
    sum(deseq2_df$log2FoldChange[deseq2_df$gene_symbol %in% consensus_genes] > 0, na.rm = TRUE)
  ),
  Downregulated = c(
    sum(deseq2_sig$log2FoldChange < 0),
    sum(limma_sig$log2FoldChange < 0),
    sum(edger_sig$log2FoldChange < 0),
    sum(deseq2_df$log2FoldChange[deseq2_df$gene_symbol %in% consensus_genes] < 0, na.rm = TRUE)
  )
)

print("Results Summary:")
print(results_summary)

# Key pathways related to cancer
cancer_pathways <- c("Cell cycle", "DNA repair", "Apoptosis", "Cell adhesion", 
                     "Immune response", "Angiogenesis", "Metabolic pathways")

if(nrow(as.data.frame(go_bp)) > 0) {
  cancer_related_go <- as.data.frame(go_bp) %>%
    filter(grepl(paste(cancer_pathways, collapse = "|"), Description, ignore.case = TRUE))
  
  cat("\nCancer-related GO terms found:", nrow(cancer_related_go), "\n")
  if(nrow(cancer_related_go) > 0) {
    print(cancer_related_go[1:5, c("Description", "pvalue", "qvalue")])
  }
}

# Export results
write.csv(consensus_genes, "consensus_DEGs.csv", row.names = FALSE)
write.csv(as.data.frame(go_bp), "GO_biological_process_enrichment.csv")
write.csv(as.data.frame(kegg_pathways), "KEGG_pathway_enrichment.csv")
write.csv(results_summary, "analysis_summary.csv")

cat("\nAnalysis completed! Key files exported:\n")
cat("- consensus_DEGs.csv: List of consensus differentially expressed genes\n")
cat("- GO_biological_process_enrichment.csv: GO enrichment results\n")
cat("- KEGG_pathway_enrichment.csv: KEGG pathway enrichment results\n")
cat("- analysis_summary.csv: Summary of all methods\n")



