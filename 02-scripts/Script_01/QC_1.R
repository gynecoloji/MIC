# ==============================================================================
# Quality Control
# ==============================================================================
# This script performs quality control of scRNA-seq data from three different samples:
# AS, AS3D, and ASPDX, followed by filtering,
# extraction of cancer cells, and visualization.
# ==============================================================================

# Load required libraries --------------------------
required_packages <- c(
  "Seurat", "dplyr", "stringr", "ggplot2", "future", "harmony", 
  "pheatmap", "RColorBrewer", "ggrepel", "SingleR", "patchwork", 
  "grid", "gridExtra", "pals", "viridis"
)

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load libraries
invisible(lapply(required_packages, library, character.only = TRUE))

# specify folders ---------------------
data_save = "../../03-data/Data_01/Processed"
data_read = "../../03-data/Data_01/Processed"
analysis_save = "../../04-analysis/Analysis_01/Ana/Topic1"
set.seed(1234)  # For reproducibility

# Define color palette --------------------------
prepare_colors <- function(n) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  sample(col_vector, n)
}

# Utility Functions --------------------------

#' Quality Control Visualization
#' 
#' @param seurat_obj Seurat object
#' @param sample_name Name of the sample for plot titles
#' @param output_dir Directory to save plots
#' 
#' @return The Seurat object with QC metrics added
plot_qc_metrics <- function(seurat_obj, sample_name, output_dir) {
  # Create output directory if it doesn't exist
  dir.create(file.path(output_dir, sample_name), showWarnings = FALSE, recursive = TRUE)
  
  # Add QC metrics
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA)/log10(seurat_obj$nCount_RNA)
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  
  df <- seurat_obj@meta.data
  df <- df %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
  df$sample <- sample_name
  
  # UMI count distribution
  p1 <- ggplot(df, aes(color=sample, x=nUMI, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  ggsave(file.path(output_dir, sample_name, paste0("QC_Cell_Density_UMI.pdf")),
         plot = p1, width = 8, height = 7)
  
  # Gene count distribution
  p2 <- df %>% 
    ggplot(aes(color=sample, x=nGene, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  ggsave(file.path(output_dir, sample_name, paste0("QC_Cell_Density_genes.pdf")),
         plot = p2, width = 8, height = 7)
  
  # Mitochondrial ratio distribution
  p3 <- df %>% 
    ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 25)
  ggsave(file.path(output_dir, sample_name, paste0("QC_Cell_Density_mitoRatio.pdf")),
         plot = p3, width = 8, height = 7)
  
  # Scatter plot of nUMI vs nGene colored by mitoRatio
  p4 <- df %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point(aes(size=log10GenesPerUMI)) + 
    scale_colour_gradient(low = "gray90", high = "black") +
    geom_smooth(method = "lm") +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 300)
  ggsave(file.path(output_dir, sample_name, paste0("QC_Cell_correlations.pdf")),
         plot = p4, width = 8, height = 7)
  
  # Gene complexity
  p5 <- df %>%
    ggplot(aes(x=log10GenesPerUMI, color=sample, fill=sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  ggsave(file.path(output_dir, sample_name, paste0("QC_Cell_complexity.pdf")),
         plot = p5, width = 8, height = 7)
  
  return(seurat_obj)
}

#' Filter cells and genes based on QC metrics
#' 
#' @param seurat_obj Seurat object with QC metrics
#' @param sample_name Name of the sample
#' @param output_dir Directory to save plots
#' 
#' @return Filtered Seurat object
filter_cells_and_genes <- function(seurat_obj, sample_name, output_dir) {
  # Cell-level filtering
  seurat_filtered <- subset(
    x = seurat_obj, 
    subset = (nCount_RNA >= 500) & 
      (nFeature_RNA >= 300) & 
      (log10GenesPerUMI > 0.80) & 
      (mitoRatio < 25)
  )
  
  # Gene-level filtering
  counts <- GetAssayData(object = seurat_filtered, slot = "counts")
  nonzero <- counts > 0
  
  # Plot gene expression distribution
  df <- data.frame(
    Gene_in_nCells = Matrix::rowSums(nonzero),
    sample = sample_name
  )
  
  p <- df %>% 
    ggplot(aes(color=sample, x=Gene_in_nCells, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + 
    geom_vline(xintercept = 10)
  ggsave(file.path(output_dir, sample_name, paste0("QC_genes_density.pdf")),
         plot = p, width = 8, height = 7)
  
  # Keep genes expressed in at least 10 cells
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  
  # Create new Seurat object with filtered data
  seurat_filtered <- CreateSeuratObject(
    filtered_counts, 
    meta.data = seurat_filtered@meta.data
  )
  
  return(seurat_filtered)
}

#' Process and cluster Seurat object
#' 
#' @param seurat_object Seurat object
#' @param npca Number of PCs to compute
#' @param ndim Number of dimensions to use for clustering
#' @param resolution Resolution parameter for clustering
#' @param var_to_regress Variables to regress out
#' 
#' @return Processed Seurat object
process_seurat <- function(seurat_object, npca=50, ndim=50, resolution=0.5, 
                          var_to_regress=c("G2M.Score","S.Score","nCount_RNA","percent.mt")) {
  
  # Add mitochondrial percentage
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # Normalize data
  seurat_object <- NormalizeData(
    seurat_object, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
  
  # Find variable features
  seurat_object <- FindVariableFeatures(
    seurat_object, 
    selection.method = "vst", 
    nfeatures = 3000
  )
  
  # Cell cycle scoring
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat_object <- CellCycleScoring(
    seurat_object, 
    s.features = s.genes, 
    g2m.features = g2m.genes, 
    set.ident = TRUE
  )
  
  # Scale data
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(
    seurat_object, 
    features = all.genes, 
    vars.to.regress = var_to_regress
  )
  
  # Run PCA
  seurat_object <- RunPCA(
    seurat_object, 
    features = VariableFeatures(object = seurat_object), 
    npcs = npca
  )
  
  # Determine significant PCs
  seurat_object <- JackStraw(seurat_object, num.replicate = 100, dims = npca)
  seurat_object <- ScoreJackStraw(seurat_object, dims = 1:npca)
  
  # Find neighbors and clusters
  seurat_object <- FindNeighbors(seurat_object, dims = 1:ndim)
  seurat_object <- FindClusters(seurat_object, resolution = resolution)
  
  # Run dimensionality reduction
  seurat_object <- RunUMAP(seurat_object, dims = 1:ndim)
  seurat_object <- RunTSNE(seurat_object, dims = 1:ndim)
  
  return(seurat_object)
}

#' Plot UMAP with labeled clusters
#' 
#' @param seurat_obj Processed Seurat object
#' @param sample_name Name of the sample
#' @param output_dir Directory to save plots
#' @param colors Color vector for clusters
#' 
#' @return NULL
plot_labeled_umap <- function(seurat_obj, sample_name, output_dir, colors) {
  x_min <- min(Embeddings(seurat_obj, "umap")[,1]) 
  y_min <- min(Embeddings(seurat_obj, "umap")[,2])
  
  # Prepare cluster label coordinates
  umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
  umap_coords$cluster <- seurat_obj@meta.data$seurat_clusters
  cluster_centers <- umap_coords %>%
    group_by(cluster) %>%
    summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
  
  # Create UMAP plot with labels
  p <- DimPlot(seurat_obj, pt.size = 1) +
    scale_color_manual(values = colors) +
    geom_label_repel(
      inherit.aes = FALSE, 
      size = 7,
      data = cluster_centers, 
      aes(x=UMAP_1, y=UMAP_2, label=cluster, fontface="bold", color="black"), 
      max.overlaps = Inf, 
      force = 8,
      max.time = 10,
      segment.size = 1,
      show.legend = FALSE
    ) +
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
      plot.title.position = "plot"
    ) +
    labs(x="UMAP 1", y="UMAP 2") +
    ggtitle(paste0("UMAP of ", sample_name)) +
    geom_segment(
      aes(x = x_min-1, y = y_min-1, xend = x_min+1, yend = y_min-1),
      arrow = arrow(length = unit(10, "pt"), type = "closed")
    ) +
    geom_segment(
      aes(x = x_min-1, y = y_min-1, xend = x_min-1, yend = y_min+1),
      arrow = arrow(length = unit(10, "pt"), type = "closed")
    ) +
    coord_fixed()
  
  ggsave(
    file.path(output_dir, sample_name, paste0("UMAP_plot.pdf")),
    plot = p, width = 8, height = 7
  )
}

#' Run SingleR annotation
#'
#' @param seurat_obj Seurat object
#' @param sample_name Name of sample
#' @param output_dir Directory to save plots
#' @param cluster_colors Colors for clusters
#'
#' @return NULL
run_singleR_annotation <- function(seurat_obj, sample_name, output_dir, cluster_colors) {
  # Get reference data
  hpca.se <- HumanPrimaryCellAtlasData()
  
  # Run SingleR
  pred.hesc <- SingleR(
    test = as.matrix(seurat_obj@assays$RNA@data), 
    ref = hpca.se,
    labels = hpca.se$label.main
  )
  
  # Prepare data for heatmap
  df <- pred.hesc$scores
  rownames(df) <- colnames(seurat_obj@assays$RNA@data)
  df_with_clusters <- as.data.frame(df)
  df_with_clusters$Cluster <- seurat_obj@meta.data$seurat_clusters
  df_with_clusters <- df_with_clusters[order(df_with_clusters$Cluster),]
  
  # Plot heatmap
  pdf(
    file = file.path(output_dir, sample_name, "SingleR_heatmap.pdf"),
    width = 12, height = 14
  )
  pheatmap::pheatmap(
    t(df_with_clusters[,-ncol(df_with_clusters)]),
    show_colnames = FALSE,
    cluster_cols = FALSE,
    annotation_col = data.frame(
      row.names = rownames(df_with_clusters), 
      Cluster = df_with_clusters$Cluster
    ),
    annotation_colors = list(Cluster = cluster_colors),
    scale = "column",
    color = viridis::viridis(20),
    fontsize = 20
  )
  dev.off()
  
  # Return table of annotations
  return(table(pred.hesc$pruned.labels))
}

# Main Function for Processing a Sample --------------------

#' Process a single sample
#'
#' @param sample_id ID of the sample
#' @param data_dir Base directory for data
#' @param output_dir Directory for outputs
#' @param is_pdx Whether sample is PDX (TRUE/FALSE)
#'
#' @return Processed Seurat object
process_sample <- function(sample_id, data_dir, output_dir, is_pdx=FALSE) {
  message(paste0("Processing sample: ", sample_id))
  
  # Data loading
  sample_path <- file.path(data_dir, sample_id, "filtered_feature_bc_matrix")
  data <- Read10X(data.dir = sample_path)
  
  if(is_pdx) {
    # Special handling for PDX sample to filter out mouse cells
    data <- as.data.frame(data)
    
    # Read gem classification file
    gem_class <- read.csv(file.path(data_dir, sample_id, "gem_classification.csv"))
    gem_class$barcode <- paste0(sample_id, "_", gem_class$barcode)
    
    # Filter for human cells only
    human_cells <- gem_class[gem_class$call %in% c("GRCh38"), "barcode"]
    data <- data[grep("^GRCh38_", rownames(data)), colnames(data) %in% human_cells]
    
    # Clean up gene names
    rownames(data) <- stringr::str_split(rownames(data), "_", n=2, simplify = TRUE)[,2]
  } else {
    # Add sample prefix to cell barcodes
    colnames(data) <- paste0(sample_id, "_", colnames(data))
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = as.matrix(data), project = sample_id)
  
  # QC and filtering
  seurat_obj <- plot_qc_metrics(seurat_obj, sample_id, output_dir)
  seurat_obj <- filter_cells_and_genes(seurat_obj, sample_id, output_dir)
  
  # Process data and generate plots
  if(sample_id == "AS3D") {
    # Full processing for AS3D with clustering
    seurat_obj <- process_seurat(
      seurat_obj,
      npca = 50,
      ndim = 50,
      resolution = 0.5,
      var_to_regress = c("G2M.Score", "S.Score", "nCount_RNA", "percent.mt")
    )
    
    # Generate cluster colors
    cluster_colors <- prepare_colors(length(levels(seurat_obj)))
    names(cluster_colors) <- levels(seurat_obj)
    
    # Plot UMAP
    plot_labeled_umap(seurat_obj, sample_id, output_dir, cluster_colors)
    
    # Plot marker genes
    p <- FeaturePlot(seurat_obj, features = c("CALB2", "EPCAM", "ACTA2", "PDPN"))
    ggsave(file.path(output_dir, sample_id, "Markers_scatterplot.pdf"),
           plot = p, width = 12, height = 12)
    
    # Cell type annotation with SingleR
    annotations <- run_singleR_annotation(seurat_obj, sample_id, output_dir, cluster_colors)
    
    # Add cell type categories
    seurat_obj@meta.data$Category <- ifelse(
      seurat_obj@meta.data$seurat_clusters == 4, 
      "Stromal_Cell", 
      "Cancer_Cell"
    )
    
    # Extract cancer cells for separate analysis
    Idents(seurat_obj) <- "Category"
    seurat_cancer_cells <- subset(seurat_obj, idents = "Cancer_Cell")
    
    # Save objects
    save(seurat_obj, file = paste0(sample_id, "_seuratObj.RData"))
    save(seurat_cancer_cells, file = paste0(sample_id, "_CC_seuratObj.RData"))
    
    return(list(full = seurat_obj, cancer_cells = seurat_cancer_cells))
  } else {
    # For AS and ASPDX, just save the filtered object
    seurat_obj@meta.data$Category <- "Cancer_Cell"
    save(seurat_obj, file = paste0(sample_id, "_seuratObj.RData"))
    return(seurat_obj)
  }
}

# Main script execution --------------------

# Create main output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Process all samples
message("Starting single-cell RNA-seq analysis pipeline")

# Process AS3D sample
as3d_results <- process_sample("AS3D", DATA_DIR, OUTPUT_DIR)

# Process AS sample 
as_results <- process_sample("AS", DATA_DIR, OUTPUT_DIR)

# Process ASPDX sample (PDX model)
aspdx_results <- process_sample("ASPDX", DATA_DIR, OUTPUT_DIR, is_pdx=TRUE)

message("Analysis completed successfully!")