#=========================================#
# Single-Cell Analysis Pipeline
# 
# This script performs several analyses on integrated scRNA-seq data:
# 1. CytoTRACE analysis for cellular trajectory inference
# 2. MAGIC imputation for gene expression visualization
# 3. Metabolic flux analysis using scFEA outputs
#=========================================#

#----------------------#
# Setup and Configuration
#----------------------#

# Python configuration (if needed)
reticulate::use_python('/Users/gynecoloji/opt/anaconda3/bin/python')

# Helper function to load required packages
load_packages <- function(pkg_list) {
  not_installed <- pkg_list[!pkg_list %in% installed.packages()[,"Package"]]
  if(length(not_installed) > 0) {
    message("Installing missing packages: ", paste(not_installed, collapse=", "))
    install.packages(not_installed)
  }
  
  invisible(lapply(pkg_list, library, character.only = TRUE))
  message("All required packages loaded successfully.")
}

# Define the required packages
required_packages <- c(
  "ggplot2", "pals", "stringr", "GSEABase", "AUCell", "RColorBrewer", 
  "ggrepel", "CytoTRACE", "Seurat", "Rmagic", "pheatmap", 
  "singleCellNet", "Matrix", "dplyr", "tidyr"
)

# Load all required packages
load_packages(required_packages)

# Set project paths
project_dir <- "../.."  # Change to your project directory
data_dir <- file.path(project_dir, "03-data/Data_01/Processed")
output_dir <- file.path(project_dir, "04-analysis/Analysis_01/Ana/Topic1/characterization")
scfea_dir <- file.path(output_dir, "metabolism/scFEA")
cytotrace_dir <- file.path(output_dir, 'Cytotrace')
hallmark_dir <- file.path(output_dir, 'HALLMARK')
markers_dir <- file.path(output_dir, 'Markers')


# Create output directories if they don't exist
dir.create(file.path(output_dir, "Cytotrace"), recursive = TRUE, showWarnings = FALSE)
dir.create(scfea_dir, recursive = TRUE, showWarnings = FALSE)


# Load the Seurat object
seurat_file <- file.path(data_dir, "Integrated_seurat_AS_AS3D_ASPDX_HVG5000.RData")
if(!file.exists(seurat_file)) {
  stop("Seurat object file not found. Please check the path: ", seurat_file)
}
load(seurat_file)

#----------------------#
# Helper Functions
#----------------------#

# Function to create consistent ggplot theme
theme_sc <- function() {
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

# Function for boxplot theme
theme_boxplot <- function() {
  theme_bw() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", arrow = arrow(length = unit(10, "pt"), type = "closed")),
      text = element_text(size = 20)
    )
}

# Function to calculate gap positions for heatmap
calculate_gaps <- function(categories) {
  # Get frequency table
  freq_table <- table(categories)
  
  # Calculate positions
  positions <- numeric(length(freq_table))
  running_sum <- 0
  
  for(i in 1:length(freq_table)) {
    running_sum <- running_sum + freq_table[i]
    positions[i] <- running_sum
  }
  
  return(positions)
}

# Function to save boxplots for gene expression
save_gene_boxplot <- function(data, gene_name, colors, output_file) {
  p <- ggplot(data, aes(x=new_type, y=!!sym(gene_name), fill=new_type)) +
    geom_boxplot() +
    scale_fill_manual(values = colors) +
    theme_boxplot() +
    xlab("") +
    ylab(gene_name)
  
  ggsave(output_file, p, width = unit(10,"inch"), height = unit(4,"inch"))
  return(p)
}

#----------------------#
# Color Preparation
#----------------------#

# Generate colors for cell types
CC_colors <- alphabet()[1:length(levels(immune.combined$new_type))]
names(CC_colors) <- paste0("CC", seq(0, length(levels(immune.combined$new_type))-1))

#----------------------#
# CytoTRACE Analysis
#----------------------#

# Prepare data for CytoTRACE
message("Running CytoTRACE analysis...")
df <- as.matrix(immune.combined@assays$integrated@data)
df <- df - min(df)  # Normalize to non-negative values
print(paste("Data range:", paste(range(df), collapse=" to ")))

# Run CytoTRACE
cytotrace_results <- CytoTRACE(df)
cell_ranks <- cytotrace_results$CytoTRACErank

# Get UMAP coordinates and add ranks
umap_coords <- as.data.frame(Embeddings(immune.combined, "umap"))
if(!identical(rownames(umap_coords), names(cell_ranks))) {
  warning("Cell order mismatch between UMAP coordinates and CytoTRACE ranks")
}
umap_coords$ranks <- cell_ranks

# Calculate axis limits for arrow visualization
x_min <- min(umap_coords$UMAP_1) 
y_min <- min(umap_coords$UMAP_2)

# Create and save UMAP plot colored by CytoTRACE ranks
p_umap <- ggplot(umap_coords, aes(x=UMAP_1, y=UMAP_2, color=ranks)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_sc() +
  labs(x="UMAP 1", y="UMAP 2") +
  ggtitle("UMAP of Cancer Cells") +
  geom_segment(aes(x = x_min-1, y = y_min-1, xend = x_min+1, yend = y_min-1),
               arrow = arrow(length = unit(10, "pt"), type = "closed"), color="black") +
  geom_segment(aes(x = x_min-1, y = y_min-1, xend = x_min-1, yend = y_min+1),
               arrow = arrow(length = unit(10, "pt"), type = "closed"), color="black") +
  coord_fixed()

ggsave(file.path(cytotrace_dir, "UMAP.pdf"), p_umap, 
       width = unit(12,"inch"), height = unit(12,"inch"))

# Add raw CytoTRACE values and cell types
umap_coords$order_values <- cytotrace_results$CytoTRACE
umap_coords$new_type <- immune.combined@meta.data$new_type

# Create rank boxplot by cell type (optional)
p_ranks <- ggplot(umap_coords, aes(x=new_type, y=ranks, fill=new_type)) +
  geom_boxplot() +
  scale_fill_manual(values = CC_colors)

#----------------------#
# MAGIC Imputation
#----------------------#

message("Running MAGIC imputation...")
# Normalize data
df <- as.matrix(immune.combined@assays$integrated@data)
df <- df - min(df)

# Run MAGIC imputation
MAGIC_data <- magic(t(df), knn = 20, t=3)
df_magic <- MAGIC_data$result
df_magic$new_type <- immune.combined$new_type

# Generate and save boxplots for key genes
gene_list <- c("MYC", "CD44", "CD24", "EPCAM", 'VIM')

for(gene in gene_list) {
  if(gene %in% colnames(df_magic)) {
    output_file <- file.path(markers_dir, paste0(gene, "_boxplot.pdf"))
    save_gene_boxplot(df_magic, gene, CC_colors, output_file)
  } else {
    warning(paste("Gene", gene, "not found in imputed data"))
  }
}

# Create ALDH aggregate boxplot
aldh_genes <- colnames(df_magic)[grep("ALDH", colnames(df_magic))]
if(length(aldh_genes) > 0) {
  tmp <- df_magic[, c(aldh_genes, "new_type")]
  tmp$ALDH <- apply(tmp[, aldh_genes], 1, mean)
  
  save_gene_boxplot(
    tmp, "ALDH", CC_colors,
    file.path(markers_dir, "ALDH_boxplot.pdf")
  )
}

#----------------------#
# scFEA Analysis - M171
#----------------------#

message("Processing scFEA M171 results...")

# Function to import and process scFEA results
process_scfea <- function(model_name, flux_file, metadata_file, celltype_data) {
  # Create output directory if it doesn't exist
  dir.create(file.path(scfea_dir, model_name), recursive = TRUE, showWarnings = FALSE)
  
  # Import flux data
  df_flux <- tryCatch({
    if(model_name == "M171") {
      df <- read.csv(flux_file, header = TRUE, row.names = 1, skip = 1)
      colnames(df) <- colnames(read.csv(flux_file, header = TRUE, row.names = 1)[1,])
      df
    } else {
      read.csv(flux_file, header = TRUE, row.names = 1)
    }
  }, error = function(e) {
    stop(paste("Error reading flux file for", model_name, ":", e$message))
  })
  
  # Import metadata
  df_meta <- read.csv(metadata_file, header = TRUE)
  
  # Check cell name consistency
  if(model_name == "M171") {
    if(!identical(rownames(df_flux), colnames(immune.combined))) {
      warning("Cell names in flux data don't match Seurat object")
    }
  }
  
  # Process flux data
  df_flux <- na.omit(t(df_flux))
  df_flux <- as.data.frame(t(df_flux))
  df_flux$new_type <- celltype_data
  df_flux <- df_flux[order(df_flux$new_type),]
  
  # Calculate average flux by cell type
  df_flux_avg <- aggregate(.~new_type, df_flux, mean)
  df_flux_avg <- df_flux_avg %>% tibble::column_to_rownames(var = "new_type")
  
  # Scale data for heatmap
  df_flux_avg_heatmap <- apply(df_flux_avg, 2, scale)
  rownames(df_flux_avg_heatmap) <- rownames(df_flux_avg)
  df_flux_avg_heatmap <- na.omit(t(df_flux_avg_heatmap))
  
  # Process annotation
  if(model_name == "M171") {
    annotation_row <- merge(df_flux_avg_heatmap, df_meta, by.x="row.names", by.y="M_id")
    class_col <- "Class"
  } else if(model_name %in% c("Lipid", "MGF")) {
    annotation_row <- merge(df_flux_avg_heatmap, df_meta, by.x="row.names", by.y="M_id")
    annotation_row$Class <- annotation_row$SL
    class_col <- "Class"
  } else {
    annotation_row <- data.frame(row.names = rownames(df_flux_avg_heatmap))
    class_col <- NULL
  }
  
  if(!is.null(class_col) && class_col %in% colnames(annotation_row)) {
    # Order by class
    annotation_row$Class <- factor(annotation_row$Class, levels = unique(annotation_row$Class))
    annotation_row <- annotation_row[order(annotation_row$Class),]
    
    # Match order in heatmap data
    df_flux_avg_heatmap <- df_flux_avg_heatmap[match(annotation_row$Row.names, rownames(df_flux_avg_heatmap)),]
    
    # Create color palette for classes
    my_colors <- alphabet2(length(unique(annotation_row$Class)))
    names(my_colors) <- unique(annotation_row$Class)
    
    # Calculate gaps for heatmap
    gaps <- calculate_gaps(annotation_row$Class)
    
    # Create heatmap with class annotation
    pdf(file.path(scfea_dir, paste0(model_name, "_modules_heatmap_avg.pdf")))
    pheatmap(
      as.matrix(df_flux_avg_heatmap),
      color = rev(brewer.pal(11,"RdBu")),
      cluster_cols = TRUE,
      cluster_rows = FALSE,
      show_colnames = TRUE,
      show_rownames = (nrow(df_flux_avg_heatmap) < 50),  # Show rownames only if not too many
      annotation_row = annotation_row[, class_col, drop=FALSE],
      annotation_colors = list(Class = my_colors),
      gaps_row = gaps,
      border_color = "black"
    )
    dev.off()
  } else {
    # Create simple heatmap without class annotation
    pdf(file.path(scfea_dir, paste0(model_name, "_modules_heatmap_avg.pdf")), 
        height = ifelse(nrow(df_flux_avg_heatmap) > 50, 20, 12), 
        width = 8)
    pheatmap(
      as.matrix(df_flux_avg_heatmap),
      color = rev(brewer.pal(11,"RdBu")),
      cluster_cols = TRUE,
      cluster_rows = FALSE,
      show_colnames = TRUE,
      show_rownames = (nrow(df_flux_avg_heatmap) < 50),
      border_color = "black"
    )
    dev.off()
  }
  
  # Return processed data
  return(list(flux = df_flux, avg_flux = df_flux_avg, annotation = annotation_row))
}

# Process M171 model
m171_results <- process_scfea(
  model_name = "M171",
  flux_file = file.path(scfea_dir, "web_output/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"),
  metadata_file = file.path(scfea_dir, "web_output/df_metaM171.csv"),
  celltype_data = immune.combined@meta.data$new_type
)

# If specific modules are of interest (M42, M58)
if(all(c("M_42", "M_58") %in% colnames(m171_results$flux))) {
  df_heatmap <- m171_results$flux[, c("M_42", "M_58", "new_type")]
  df_heatmap <- df_heatmap[order(df_heatmap$new_type),]
  
  # Calculate ALDH associated module scores
  df_heatmap_scale <- data.frame(
    ALDH_associated_Modules = apply(df_heatmap[, c("M_42", "M_58")], 1, sum),
    new_type = df_heatmap$new_type
  )
  
  # Create boxplot
  p_aldh_modules <- ggplot(df_heatmap_scale, aes(x=new_type, y=ALDH_associated_Modules, fill=new_type)) +
    geom_boxplot() +
    scale_fill_manual(values = CC_colors) +
    theme_boxplot() +
    xlab("") + 
    ylab("ALDH associated Modules")
  
  ggsave(file.path(output_dir, "ALDH_associated_Modules_boxplot.pdf"), 
         p_aldh_modules, width = unit(10,"inch"), height = unit(4,"inch"))
}

#----------------------#
# Process other scFEA models
#----------------------#

# Process Lipid model
if(file.exists(file.path(scfea_dir, "lipid/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"))) {
  process_scfea(
    model_name = "Lipid",
    flux_file = file.path(scfea_dir, "lipid/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"),
    metadata_file = file.path(scfea_dir, "lipid/scFEA.Lipid-metabolism.human.moduleinfo.csv"),
    celltype_data = immune.combined@meta.data$new_type
  )
}

# Process MGF model
if(file.exists(file.path(scfea_dir, "MGF/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"))) {
  process_scfea(
    model_name = "MGF",
    flux_file = file.path(scfea_dir, "MGF/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"),
    metadata_file = file.path(scfea_dir, "MGF/scFEA.MGF.human.moduleinfo.csv"),
    celltype_data = immune.combined@meta.data$new_type
  )
}

# Process KEGG model
if(file.exists(file.path(scfea_dir, "KEGG/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"))) {
  process_scfea(
    model_name = "KEGG",
    flux_file = file.path(scfea_dir, "KEGG/AS_AS3D_ASPDX_integrated_normalized_data_flux.csv"),
    metadata_file = file.path(scfea_dir, "KEGG/scFEA.KEGG.human.moduleinfo.csv"),
    celltype_data = immune.combined@meta.data$new_type
  )
}

message("Analysis complete! All results saved to output directories.")