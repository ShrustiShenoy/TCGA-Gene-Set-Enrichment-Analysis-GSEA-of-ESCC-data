# TCGA-Gene-Set-Enrichment-Analysis-GSEA-of-ESCC-data
This Python script performs Gene Set Enrichment Analysis (GSEA) across different cancer stages (Stage I–IV) using mutation data in MAF format. It extracts genes with impactful variants and identifies enriched Gene Ontology (GO) terms using the Enrichr database.

# Directory Structure
.

├── tcga_gsea3.py             # Main script

├── grade_generalised/        # Folder containing MAF files organized by stage and case_id

│   ├── StageI/

│   ├── StageII/

│   ├── StageIII/

│   └── StageIV/

└── gsea_result_all_genes/    # Output results (plots, Excel, heatmap)

# Working
1.Loads .maf files for each stage and filters impactful mutations.

2.Extracts unique genes per stage.

3.Performs GSEA using gseapy's Enrichr interface for:

GO Biological Process (BP)

GO Molecular Function (MF)

GO Cellular Component (CC)

# Output

Top 10 enriched GO terms (based on adjusted p-value)

Barplots for each stage & ontology

A cross-stage heatmap based on Combined Score

Excel file summarizing results for all stages
