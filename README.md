# 🧬 Colorectal Cancer Differential Gene Expression & Network Analysis

This project identifies robust differentially expressed genes (DEGs) in colorectal cancer using DESeq2, limma, and edgeR. It performs functional enrichment and protein–protein interaction (PPI) analysis to identify key genes and biological processes.

## 📌 Highlights

- 🎯 **Consensus DEGs** across DESeq2, limma, and edgeR: **3,566 genes**
- 🔍 **GO terms enriched**: Immune-related processes
- 🔗 **Hub genes**: MYC, VEGFA, CDK1, MMP9, IL6, FN1, CD44, etc.
- 🧠 **Biological focus**: Tumor microenvironment, inflammation, angiogenesis

## 📂 Project Structure

- `data/`: Raw and processed input data
- `results/`: Exported CSVs (DEGs, enrichment, summary)
- `plots/`: Publication-ready visualizations
- `scripts/`: All analysis scripts in R

## 📊 Results Summary

| Method   | DEGs | Upregulated | Downregulated |
|----------|------|-------------|---------------|
| DESeq2   | 5803 | 3246        | 2557          |
| limma    | 8047 | 5048        | 2999          |
| edgeR    | 3659 | 2029        | 1630          |
| Consensus| 3566 | 1982        | 1584          |

## 🔬 Enrichment Highlights

- **GO:** Humoral immune response, adaptive immune system, cytokine activity
- **KEGG:** NF-κB signaling, cytokine–cytokine receptor interaction

## 🔗 PPI Network

- Nodes: 2667  
- Edges: 69,788  
- Top hub genes: `MYC`, `VEGFA`, `CDK1`, `MMP9`, `IL6`, `FN1`, `CD44`
