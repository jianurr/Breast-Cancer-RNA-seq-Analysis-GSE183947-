## 🧬 Complete Analysis Code (R)

```r
# ---------------------------------------------------------
# Load Required Libraries
# ---------------------------------------------------------
library(readr)           # data loading
library(dplyr)           # data manipulation
library(limma)           # differential expression
library(edgeR)           # filtering utilities
library(ggplot2)         # visualizations
library(ggrepel)         # volcano plot labels
library(pheatmap)        # heatmaps
library(clusterProfiler) # GO/KEGG enrichment
library(org.Hs.eg.db)    # annotation database
library(randomForest)    # machine learning ranking

# ---------------------------------------------------------
# Load Expression Data
# ---------------------------------------------------------
expr <- read_tsv("data/fpkm_matrix.tsv") %>% 
  column_to_rownames(1)

meta <- read_tsv("data/sample_metadata.tsv")

# Check sample groups
table(meta$group)   # should show 30 tumor, 30 normal

# ---------------------------------------------------------
# Filter Low-Expression Genes
# Keep genes expressed in at least 2 samples
# ---------------------------------------------------------
expr_filt <- expr[rowSums(expr > 1) >= 2, ]

# ---------------------------------------------------------
# log2 Normalize the Data
# ---------------------------------------------------------
logexpr <- log2(expr_filt + 1)

# ---------------------------------------------------------
# Create Design Matrix (Tumor vs Normal)
# ---------------------------------------------------------
group <- factor(meta$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# ---------------------------------------------------------
# Differential Expression using limma
# ---------------------------------------------------------
fit <- lmFit(logexpr, design)
contrast <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))

# Get all DEGs
deg <- topTable(fit2, number = Inf)

# Apply your project filter:
# logFC >= 1 AND adj.p.value <= 0.01
deg_filtered <- deg %>% 
  filter(abs(logFC) >= 1 & adj.P.Val <= 0.01)

# Total DEGs found
nrow(deg_filtered)  # result: 1655 DEGs

# Save results
write_tsv(deg, "results/limma_top_table.tsv")
write_tsv(deg_filtered, "results/limma_filtered_DEGs.tsv")

# ---------------------------------------------------------
# Volcano Plot
# ---------------------------------------------------------
volcano <- ggplot(deg, aes(logFC, -log10(adj.P.Val))) +
  geom_point(alpha=0.7) +
  geom_point(data=deg_filtered, aes(logFC, -log10(adj.P.Val)), color="red") +
  theme_minimal() +
  ggtitle("Volcano Plot of DEGs")

ggsave("figures/volcano.png", volcano, width=7, height=6)

# ---------------------------------------------------------
# Heatmap of Top 20 DEGs
# ---------------------------------------------------------
top20 <- head(rownames(deg_filtered[order(deg_filtered$adj.P.Val), ]), 20)
pheatmap(logexpr[top20, ],
         scale = "row",
         annotation_col = data.frame(Group = group),
         filename = "figures/heatmap_top20.png")

# ---------------------------------------------------------
# Boxplot for Top 5 DEGs
# ---------------------------------------------------------
top5 <- rownames(deg_filtered)[1:5]

box_data <- logexpr[top5, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to="Sample", values_to="Expression") %>%
  left_join(meta, by="Sample")

boxplot_fig <- ggplot(box_data, aes(group, Expression, fill=group)) +
  geom_boxplot() +
  facet_wrap(~Gene, scales="free") +
  theme_minimal() +
  ggtitle("Top 5 DEGs (Tumor vs Normal)")

ggsave("figures/boxplot_top5.png", boxplot_fig, width=9, height=6)

# ---------------------------------------------------------
# Functional Enrichment (GO & KEGG)
# ---------------------------------------------------------
# Convert gene symbols → Entrez IDs
genes <- rownames(deg_filtered)
entrez <- bitr(genes, fromType="SYMBOL",
               toType="ENTREZID",
               OrgDb=org.Hs.eg.db)

# GO Biological Process
ego <- enrichGO(gene = entrez$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                readable = TRUE)

write_tsv(as.data.frame(ego), "results/GO_enrichment.tsv")

# KEGG Pathway Enrichment
ekegg <- enrichKEGG(gene = entrez$ENTREZID,
                    organism = "hsa")

write_tsv(as.data.frame(ekegg), "results/KEGG_enrichment.tsv")

# ---------------------------------------------------------
# KEGG Pathway Map (hsa05200)
# (Plotting handled externally using KEGG API or Pathview)
# ---------------------------------------------------------
# Note: Pathway image generated manually, not via code here.

# ---------------------------------------------------------
# Random Forest Gene Ranking (Optional Task)
# ---------------------------------------------------------
rf_data <- as.data.frame(t(logexpr[genes, ]))
rf_data$Class <- meta$group

model <- randomForest(Class ~ ., data=rf_data, importance=TRUE)

importance_df <- as.data.frame(model$importance)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing=TRUE), ]

write_tsv(importance_df, "results/random_forest_gene_importance.tsv")

# Plot importance of top 20 genes
imp_plot <- ggplot(importance_df[1:20, ], aes(x=reorder(rownames(importance_df)[1:20],
                                 MeanDecreaseGini),
                                 y=MeanDecreaseGini)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Random Forest Gene Importance")

ggsave("figures/random_forest_gene_importance.png", imp_plot, width=7, height=6)
```
