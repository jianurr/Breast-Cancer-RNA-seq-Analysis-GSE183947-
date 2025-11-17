# Breast Cancer RNA-seq Differential Expression and Functional Analysis

This project presents a complete in-silico RNA-seq differential expression and pathway analysis using the publicly available **GSE183947** breast cancer dataset. The objective of the Weblem-3 assignment is to identify differentially expressed genes (DEGs) between **breast tumor** and **normal** tissues and investigate their involvement in relevant biological functions and pathways.

The analysis is performed directly on an **FPKM expression matrix**, and therefore raw read mapping and quantification steps were skipped.

---

##  **Objective**
To perform differential expression analysis and functional enrichment on breast cancer RNA-seq data, and to explore gene involvement in biological processes, disease associations, and cancer-related pathways.

---

##  **Dataset Information**
- **Source:** GEO (GSE183947)  
- **Total Samples:** 60  
  - **30 Breast Tumor samples**  
  - **30 Normal tissue samples**  
- **Data Type:** FPKM (Fragments Per Kilobase of transcript per Million mapped reads)  
- **Genes Loaded:** 20,246  

FPKM values were used directly for downstream normalization and modeling.

---

##  **Analysis Workflow**

### **1. Data Loading**
- Imported the FPKM gene expression matrix and sample metadata.  
- Verified sample distribution: 30 tumors + 30 normals.

### **2. Preprocessing & Filtering**
- Removed low-expression genes to reduce noise.  
- Retained genes expressed in at least **2 samples**, ensuring sufficient biological signal.

---

##  **Differential Expression Analysis (limma)**
- Created a **Tumor vs Normal** design matrix.  
- Applied **linear modeling** and **empirical Bayes moderation**.  
- Applied significance thresholds:
  - **logFC ≥ 1**
  - **adj.p-value ≤ 0.01**
- **Total DEGs identified:** **1655 significant genes**

These genes represent those most strongly dysregulated in breast cancer tissue.

---

##  **Visualizations**

### **🔹 Volcano Plot**
Shows significantly up- and down-regulated genes.  
Significant DEGs appear clearly highlighted.

### **🔹 Heatmap (Top 20 DEGs)**
- Displays expression of top 20 DEGs across all 60 samples.  
- Shows clear separation and clustering of tumor vs normal groups.

### **🔹 Boxplot (Top 5 DEGs)**
- Compares expression patterns of the top 5 DEGs between tumor and normal samples.  
- Highlights strong differential expression.

---

##  **Functional Enrichment Analysis**

Performed to understand the biological roles of DEGs:

### **🔹 GO Enrichment (Biological Process)**
Revealed processes linked to:
- Cancer progression  
- Immune response  
- Cell-cycle alterations  
- Tissue remodeling  

### **🔹 KEGG Pathway Enrichment**
Identified pathways such as:
- Pathways in Cancer  
- Immune and inflammatory pathways  
- Viral oncogenesis  
- Metabolic dysregulation  

---

## 🗺 **KEGG Pathway Map – hsa05200 (“Pathways in Cancer”)**
A customized KEGG pathway map was generated:
- **Orange:** Up-regulated genes in tumor  
- **Blue:** Down-regulated genes  
- **Grey/White:** Unchanged or unmeasured  

This visualization highlights the dysregulated signaling modules in breast cancer.

---

##  **Machine Learning (Optional Task)**
A **Random Forest model** was trained using the FPKM values of top DEGs.

### Output:
- A ranked list of genes based on importance in classifying **Tumor vs Normal**.  
- A bar plot showing gene importance scores.

---

##  **Generated Outputs**
The analysis produced the following files:

### **Results**
- `limma_top_table.tsv` — complete DEG list  
- `annotated_DEGs.tsv` — DEG table with gene annotations  

### **Plots**
- `volcano.png`  
- `heatmap_top20.png`  
- `boxplot_top5.png`  
- `GO_enrichment.png`  
- `KEGG_enrichment.png`  
- `kegg_pathway_hsa05200.png`  
- `rf_gene_importance.png`

---

##  **How to Run This Project**

### 1. Place Data in `/data`
- `fpkm_matrix.tsv`  
- `sample_metadata.tsv`

### 2. Install Required Packages
```r
source("scripts/install_packages.R")
```

### 3. Run Notebook
Open:  
`notebook/weblem_analysis.Rmd`  
Run all steps to reproduce results.

---

## 👤 **Author**
**Khandakar Jianur Islam**

---

## 📄 **License**
MIT License
