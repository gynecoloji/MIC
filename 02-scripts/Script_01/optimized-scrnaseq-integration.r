# ==============================================================================
# Cancer Cell Integration and Analysis Pipeline
# ==============================================================================
# This script performs integration of scRNA-seq data from three different samples:
# AS, AS3D, and ASPDX cancer cells, followed by dimensionality reduction,
# clustering, and visualization.
# ==============================================================================

# Load required libraries --------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(SingleR)
  library(patchwork)
  library(grid)
  library(gridExtra)
  library(pals)
  library(harmony)
  library(future)
  library(msigdbr)
  library(tidyr)
})

# Set up future for parallel processing
plan(multisession, workers = 4)

# specify folders ---------------------
data_save = "../../03-data/Data_01/Processed/"
data_read = "../../03-data/Data_01/Processed/"
analysis_save = "../../04-analysis/Analysis_01/Ana/Topic1/"

# Data Loading -------------------------------------------------------------
# Load pre-processed Seurat objects
load_and_rename <- function(file_path, new_name) {
  load(file_path)
  if (exists("seurat")) {
    assign(new_name, seurat, envir = .GlobalEnv)
  } else if (exists("seurat_AS3D_CC")) {
    assign(new_name, seurat_AS3D_CC, envir = .GlobalEnv)
  }
}

load_and_rename(str_c(data_read,"ASPDX_CC_seuratObj.RData"), "seurat_ASPDX")
load_and_rename(str_c(data_read,"AS3D_CC_seuratObj.RData"), "seurat_AS3D")
load_and_rename(str_c(data_read,"AS_seuratObj.RData"), "seurat_AS")

# Clean up environment
rm(list = c("seurat", "seurat_AS3D_CC"), envir = .GlobalEnv, inherits = FALSE)

# Integration of Samples ---------------------------------------------------
# Create list of Seurat objects
seurat.list <- list(AS = seurat_AS, AS3D = seurat_AS3D, ASPDX = seurat_ASPDX)

# Normalize and find variable features for each dataset
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 5000)

# Find integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                        anchor.features = features)

# Integrate data
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# Dimensionality Reduction and Clustering ----------------------------------
# Run standard workflow for dimensionality reduction and clustering
immune.combined <- immune.combined %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(n.components = 2, reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# Label clusters with consistent naming convention
immune.combined$new_type <- paste0("CC", immune.combined$seurat_clusters)
cluster_levels <- paste0("CC", sort(as.numeric(unique(immune.combined$seurat_clusters))))
immune.combined@meta.data$new_type <- factor(immune.combined@meta.data$new_type, levels = cluster_levels)
Idents(immune.combined) <- "new_type"

# Cell Cycle Scoring ------------------------------------------------------
# Score cells for cell cycle phases
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, 
                                   g2m.features = g2m.genes, set.ident = TRUE)

# Set up colors and coordinates for visualization -------------------------
# Define color palettes
n_clusters <- length(unique(immune.combined$new_type))
CC_colors <- alphabet()[1:n_clusters]
names(CC_colors) <- sort(unique(immune.combined$new_type))

CC_colors_TS <- col_vector[71:73]
names(CC_colors_TS) <- c("AS", "AS3D", "ASPDX")

CC_colors_Phase <- col_vector[68:70]
names(CC_colors_Phase) <- c("G1", "S", "G2M")

# Set sample and phase as factors with specific levels
immune.combined@meta.data$Sample <- factor(immune.combined@meta.data$orig.ident, 
                                         levels = names(CC_colors_TS))
immune.combined@meta.data$Phase <- factor(immune.combined@meta.data$Phase, 
                                        levels = names(CC_colors_Phase))

# Calculate coordinates for text labels and arrows
x_min <- min(Embeddings(immune.combined, "umap")[,1])
y_min <- min(Embeddings(immune.combined, "umap")[,2])

# Calculate cluster centers for labels (optional)
CC_text <- as.data.frame(Embeddings(immune.combined, "umap"))
CC_text$new_type <- immune.combined@meta.data$new_type
CC_text <- aggregate(.~new_type, CC_text, mean)

# Add blank row for proper formatting
tmp <- data.frame(
  new_type = "",
  UMAP_1 = NA,
  UMAP_2 = NA
)
CC_text <- rbind(tmp, CC_text)

# Visualization Functions -------------------------------------------------
# Create a common theme for UMAP plots
umap_theme <- function() {
  theme(
    text = element_text(size = 20),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(hjust = 0.1, size = 10),
    axis.title.y = element_text(hjust = 0.1, size = 10),
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank()
  )
}

# Function to add arrows to UMAP
add_umap_arrows <- function(plot) {
  plot + 
    geom_segment(aes(x = x_min-1, y = y_min-1, xend = x_min+1, yend = y_min-1),
                 arrow = arrow(length = unit(10, "pt"), type = "closed"), color="black") +
    geom_segment(aes(x = x_min-1, y = y_min-1, xend = x_min-1, yend = y_min+1),
                 arrow = arrow(length = unit(10, "pt"), type = "closed"), color="black") +
    coord_fixed()
}

# Function for barplot theme
barplot_theme <- function() {
  theme_bw() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", 
                              arrow = arrow(length = unit(10, "pt"), type = "closed")),
      text = element_text(size = 20)
    )
}

# Visualizations ---------------------------------------------------------
# 1. UMAP by cluster
p1 <- ggplot(data = FetchData(immune.combined, c("new_type","UMAP_1","UMAP_2")), 
             aes(x = UMAP_1, y = UMAP_2, color = new_type)) +
  geom_point(alpha = .8, size = 3) +
  scale_color_manual(values = CC_colors) +
  geom_label_repel(
    inherit.aes = FALSE, 
    size = 7,
    data = CC_text, 
    aes(x = UMAP_1, y = UMAP_2, label = new_type, fontface = "bold", color = "black"),
    max.overlaps = Inf,
    force = 50,
    max.time = 100,
    segment.size = 1,
    show.legend = FALSE
  ) +
  umap_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP of Cancer Cells")

p1 <- add_umap_arrows(p1)
ggsave(str_c(analysis_save, "UMAP_of_CC.pdf"), p1, width = 12, height = 12)

# 2. UMAP by sample
p2 <- ggplot(data = FetchData(immune.combined, c("Sample","UMAP_1","UMAP_2")), 
             aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
  geom_point(alpha = .8, size = 3) +
  scale_color_manual(values = CC_colors_TS) +
  umap_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP of Cancer Cells")

p2 <- add_umap_arrows(p2)
ggsave(str_c(analysis_save,"UMAP_of_CC_in_Samples.pdf"), p2, width = 8, height = 6)

# 3. UMAP by cell cycle phase
p3 <- ggplot(data = FetchData(immune.combined, c("Phase","UMAP_1","UMAP_2")), 
             aes(x = UMAP_1, y = UMAP_2, color = Phase)) +
  geom_point(alpha = .8, size = 3) +
  scale_color_manual(values = CC_colors_Phase) +
  umap_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP of Cancer Cells")

p3 <- add_umap_arrows(p3)
ggsave(str_c(analysis_save,"UMAP_of_CC_Cycle.pdf"), p3, width = 8, height = 6)

# 4. Line plot of cluster percentages across samples
cluster_percentage_by_sample <- immune.combined@meta.data %>%
  group_by(Sample, new_type) %>%
  summarize(Count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count)) %>%
  ungroup()

p4 <- ggplot(cluster_percentage_by_sample, aes(x = Sample, y = Percentage, group = new_type)) + 
  geom_line(aes(color = new_type), linewidth = 2) +
  geom_point(aes(color = new_type), size = 4) +
  scale_color_manual(values = CC_colors) + 
  theme_bw() + 
  xlab("Sample") + 
  ylab("Percentage") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", 
                            arrow = arrow(length = unit(10, "pt"), type = "closed")),
    text = element_text(size = 10)
  ) +
  facet_wrap(~new_type, nrow = 3)

ggsave(str_c(analysis_save,"lineplot_of_CC_among_Samples.pdf"), p4, width = 16, height = 6)

# 5. Barplot of clusters by sample
p5 <- ggplot(immune.combined@meta.data, aes(fill = new_type, x = Sample)) + 
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = CC_colors) + 
  barplot_theme() +
  xlab("Sample") + 
  ylab("Percentage")

ggsave(str_c(analysis_save,"Barplot_of_CC_Samples.pdf"), p5, width = 5, height = 6)

# 6. Barplot of cell cycle by cluster
p6 <- ggplot(immune.combined@meta.data, aes(fill = Phase, x = new_type)) + 
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = CC_colors_Phase) + 
  barplot_theme() +
  xlab("Cell Clusters") + 
  ylab("Percentage")

ggsave(str_c(analysis_save,"Barplot_of_CC_Cell_Cycle.pdf"), p6, width = 10, height = 4)

# 7. Barplot of samples by cluster
p7 <- ggplot(immune.combined@meta.data, aes(fill = Sample, x = new_type)) + 
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = CC_colors_TS) + 
  barplot_theme() +
  xlab("Cell Clusters") + 
  ylab("Percentage")

ggsave(str_c(analysis_save,"Barplot_of_Sample_in_CC.pdf"), p7, width = 10, height = 4)

# Save integrated Seurat object and metadata ------------------------------
save(immune.combined, file = str_c(data_save,"Integrated_seurat_AS_AS3D_ASPDX_HVG5000.RData"))

# Generate cluster distribution matrix
cluster_distribution <- table(immune.combined$new_type, immune.combined$Sample) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  column_to_rownames(var = "Var1")

# Export results for other analyses
write.csv(t(cluster_distribution), "Distribution_CC_in_Samples.csv")

# Create heatmap of cluster distribution
pdf(str_c(analysis_save,"Distribution_CC_in_Samples.pdf"), width = 12, height = 8)
pheatmap(
  t(cluster_distribution),
  color = rev(brewer.pal(11, "RdBu")),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  display_numbers = t(cluster_distribution),
  fontsize_number = 20,
  fontsize = 20
)
dev.off()

# Export data for RNA velocity analysis (optional)
dir.create(str_c(data_save,"Velocyto"))
write.csv(immune.combined@meta.data, str_c(data_save,"Velocyto/Metadata.csv"))
write.csv(Embeddings(immune.combined, 'umap'), str_c(data_save,"Velocyto/UMAP_embedding.csv"))
write.csv(VariableFeatures(immune.combined),str_c(data_save, "Velocyto/HVG.csv"))
write.csv(CC_colors, str_c(data_save,"Velocyto/CC_colors.csv"))

# Create Multi-panel Figure -----------------------------------------------
# Use smaller text sizes for multi-panel display
f1 <- DimPlot(immune.combined, group.by = "Sample", pt.size = 1) +
  scale_color_manual(values = CC_colors_TS) +
  theme(text = element_text(size = 10)) +
  umap_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP by Sample")

f1 <- add_umap_arrows(f1)

f2 <- DimPlot(immune.combined, pt.size = 1) +
  scale_color_manual(values = CC_colors) +
  geom_label_repel(
    inherit.aes = FALSE, 
    size = 3,
    data = CC_text, 
    aes(x = UMAP_1, y = UMAP_2, label = new_type, fontface = "bold", color = "black"),
    max.overlaps = Inf, 
    force = 50,
    max.time = 100,
    segment.size = 1,
    show.legend = FALSE
  ) +
  theme(text = element_text(size = 10)) +
  umap_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP by Cluster")

f2 <- add_umap_arrows(f2)

f3 <- ggplot(immune.combined@meta.data, aes(fill = new_type, x = Sample)) + 
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = CC_colors) + 
  theme_bw() +
  xlab("Sample") +
  ylab("Percentage") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", 
                             arrow = arrow(length = unit(10, "pt"), type = "closed")),
    text = element_text(size = 10)
  )

f4 <- ggplot(immune.combined@meta.data, aes(fill = Sample, x = new_type)) + 
  geom_bar(position = "fill", stat = "count") +
  scale_fill_manual(values = CC_colors_TS) + 
  theme_bw() +
  xlab("Cluster") +
  ylab("Percentage") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", 
                             arrow = arrow(length = unit(10, "pt"), type = "closed")),
    text = element_text(size = 10)
  )

# Create layout for multi-panel figure
layout <- rbind(
  c(1, 1, 1, 1, 3, 3, 3),
  c(1, 1, 1, 1, 3, 3, 3),
  c(1, 1, 1, 1, 3, 3, 3),
  c(2, 2, 2, 2, 4, 4, 4),
  c(2, 2, 2, 2, 4, 4, 4),
  c(2, 2, 2, 2, 4, 4, 4)
)

# Save multi-panel figure
pdf("Summary_Panel_Figure.pdf", width = 14, height = 10)
grid.arrange(f1, f2, f3, f4, layout_matrix = layout)
dev.off()
