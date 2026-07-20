# Clear workspace
rm(list=ls(all=TRUE))

# # Install BiocManager if not already installed and then use it to install Mfuzz
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("Mfuzz")

# Load required libraries
invisible(suppressWarnings(suppressPackageStartupMessages(lapply(c(
  "DESeq2", "ggplot2", "ggrepel", "dplyr", "tidyr", "reshape2", "pheatmap", "matrixStats",
  "Mfuzz", "optparse", "cowplot", "argparse", "ggplotify", "pracma", "umap", "Rtsne"
), library, character.only = TRUE))))

set.seed(2025)  # Fix randomness

# Define command-line arguments
option_list <- list(
  make_option(c("-d", "--design"), type="character", help="Path to design matrix file"),
  make_option(c("-c", "--contrast"), type="character", help="Path to contrast file"),
  make_option(c("-i", "--input"), type="character", help="Directory containing RNA-seq count files"),
  make_option(c("-p", "--prefix"), type="character", help="Output prefix"),
  make_option(c("-a", "--alpha"), type="double", default=0.05, help="Padj cutoff for DEG selection"),
  make_option(c("-l", "--log2fc"), type="double", default=0.263, help="Log2FC cutoff for DEG selection"),
  make_option(c("-m", "--mode"), type="character", default="any", help="Mode: 'all' or 'any'")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$design) || is.null(opt$contrast) || is.null(opt$input) || is.null(opt$prefix)) {
  print_help(opt_parser)
  stop("All arguments must be provided.", call.=FALSE)
}

# Read design matrix and contrast file
design_matrix <- read.table(opt$design, header=TRUE, sep="\t")
design_matrix$Genotype <- as.factor(design_matrix$Genotype)
contrasts <- read.table(opt$contrast, header=FALSE, sep="\t")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=design_matrix, directory=opt$input, design=~Genotype)
dds <- estimateSizeFactors(ddsHTSeq)

# Apply row mean cutoff for filtering low-expressed genes
dat <- counts(dds, normalized=TRUE)
# Ensure no negative counts in the data
dat[dat < 0] <- 0  # Replace negative values with zero
idx <- rowMeans(dat) > 50
dds <- dds[idx,]
dat <- dat[idx,]

dds <- DESeq(dds, betaPrior=FALSE)
write.table(dat, file=paste0(opt$prefix,"_NormalizedCounts.txt"), sep="\t", row.names=TRUE, quote=FALSE)

# Select DEGs for Clustering
deg_list <- list()
for (i in 1:nrow(contrasts)) {
  res <- results(dds, contrast=c("Genotype", contrasts[i,2], contrasts[i,3]))
  res$FC <- 2^res$log2FoldChange
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  write.table(resdata, file=paste0(opt$prefix,"_", contrasts[i,1], "_DEG.txt"), sep="\t", row.names=FALSE, quote=FALSE)

  deg_list[[contrasts[i,1]]] <- resdata %>% 
    filter(padj <= opt$alpha & abs(log2FoldChange) >= opt$log2fc) %>%
    select(Row.names)
}

if (opt$mode == "all") {
  deg_genes <- Reduce(intersect, deg_list)
} else {
  deg_genes <- unique(unlist(deg_list))
}

# Extract DEG data and log-transform counts
deg_data <- log1p(dat[deg_genes, ])

# Convert to ExpressionSet
eset <- new("ExpressionSet", exprs=as.matrix(deg_data))

# Standardize expression data
eset <- standardise(eset)

# Assign time points based on sample order
time_points <- rep(1:6, each=3)  # Assuming 5 time points, 4 replicates each
pData(eset) <- data.frame(Time = time_points)

c_fixed <- 7

# Estimate optimal fuzziness parameter (m)
m_optimal <- mestimate(eset)

# Compute Dmin for a range of clusters (2 to 15)
c_range <- 2:15
dmin_values <- Dmin(eset, m=m_optimal, crange=c_range, repeats=3)


# Method 1: Relative Change Method (<10%) -- work in progress
#rel_change <- diff(dmin_values) / dmin_values[-length(dmin_values)]
#c_rel_change <- c_range[which(rel_change < 0.05)[1] + 1]

# Method 2: Slope Method
slope <- diff(dmin_values) / diff(c_range)

# Compute slope differences (first derivative)
slope_diff <- diff(slope)

# Find all indices where the slope increases (slope[i] > slope[i-1])
increase_indices <- which(slope_diff > 0) + 1  # Adjust for diff() shifting index

# Identify the last index before the first decrease after a continuous increase
stabilization_index <- NA
last_increase_index <- NA

for (i in 2:(length(slope) - 1)) {  # Start from index 2 since index 1 has no previous value
  if (slope[i] > slope[i - 1]) {  # Detect increasing sequence
    last_increase_index <- i  # Track last increasing index
  } else if (!is.na(last_increase_index)) {  # First dip detected after increase
    stabilization_index <- last_increase_index  # Return last increasing index before dip
    break
  }
}

#print(slope) 
#print(slope_diff)
#print(stabilization_index)

# Handling the case where slope never drops
if (is.na(stabilization_index)) {
  # Find the first index where the slope change is minimal (stabilization effect)
  min_slope_change_index <- which(slope_diff < 0.005)[1] + 1  # First point where slope stops increasing significantly
  
  if (!is.na(min_slope_change_index)) {
    stabilization_index <- min_slope_change_index
  } else {
    stabilization_index <- length(slope)  # Default to the last index if no stabilization is found
  }
}

# Get corresponding number of clusters
# Select best cluster numbers
#c_optimal <- unique(c(c_rel_change, c_slope))

c_optimal <- stabilization_index + 1  # Adjust to match c_range indexing

#print(c_optimal)


# Run Mfuzz clustering
#mfuzz_result <- mfuzz(eset, c=max(c_optimal), m=m_optimal)
#mfuzz_result <- mfuzz(eset, c=c_optimal, m=m_optimal,iter.max = 100000000)
mfuzz_result <- mfuzz(eset, c=c_fixed, m=m_optimal, iter.max = 100000000)


# Save cluster membership files
for (i in 1:max(c_optimal)) {
  cluster_genes <- rownames(deg_data)[mfuzz_result$cluster == i]
  write.table(cluster_genes, file=paste0(opt$prefix, "_Cluster_BeforeFiltering", i, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
}

# Save Clustering Quality Metrics
quality_metrics <- data.frame(
  Optimal_C = c_optimal,
  Optimal_m = m_optimal
)
write.table(quality_metrics, file=paste0(opt$prefix, "_ClusteringQualityMetrics.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Generate and save elxw plot for C selection


png(paste0(opt$prefix, "_Dmin_ElbowPlot.png"), width=2400, height=1800, res=300)
plot(c_range, dmin_values, type="b", pch=19, col="blue", xlab="Number of Clusters (C)", ylab="Dmin Value",
     main="Optimal Cluster Number Selection")
abline(v=c_optimal, col="red", lty=2)
dev.off()

# Extract membership matrix
membership <- mfuzz_result$membership
max_membership <- apply(membership, 1, max)
second_best <- apply(membership, 1, function(x) sort(x, decreasing=TRUE)[2])
diff12 <- max_membership - second_best

# Identify genes for removal
bottom_5_cutoff <- quantile(max_membership, 0.05)
low_confidence <- which(max_membership <= bottom_5_cutoff)
ambiguous <- which(diff12 < 0.10)
remove_idx <- union(low_confidence, ambiguous)
keep_idx <- setdiff(seq_len(nrow(membership)), remove_idx)
filtered_clusters <- mfuzz_result$cluster
filtered_clusters[remove_idx] <- 0

# Save cluster membership after filtering
for (i in unique(filtered_clusters)) {
  cluster_genes <- rownames(deg_data)[filtered_clusters == i]
  write.table(cluster_genes, file=paste0(opt$prefix, "_Cluster", i, "_AfterFiltering.txt"), sep="\t", row.names=FALSE, quote=FALSE)
}

# Ensure annotation_col is properly defined
annotation_col <- data.frame(Genotype = design_matrix$Genotype)
rownames(annotation_col) <- colnames(deg_data)

# Define sample mapping
sample_mapping <- data.frame(Sample = colnames(deg_data), Genotype = design_matrix$Genotype)



# Function to generate cluster heatmaps and trend plots
generate_cluster_plots <- function(cluster_data, cluster_labels, prefix_suffix) {
  for (clust in unique(cluster_labels)) {
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    gene_count <- length(cluster_genes)
    if (gene_count == 0) next  # Skip empty clusters

    # Ensure annotation follows the same order
    annotation_col <- data.frame(Genotype = factor(design_matrix$Genotype, levels=unique(design_matrix$Genotype)))
    rownames(annotation_col) <- colnames(cluster_data)

    # Generate heatmap
    png(paste0(opt$prefix, "_Cluster", clust, "_", prefix_suffix, "_Heatmap.png"), width=2400, height=1800, res=300)
    pheatmap(cluster_data[cluster_genes, ], scale="row", cluster_rows=FALSE, cluster_cols=FALSE,
             show_rownames=FALSE, annotation_col=annotation_col, annotation_names_col=FALSE,
             main=paste("Cluster", clust, "-", gene_count, "Genes"))
    dev.off()

    # Generate trend plots
    mean_expression <- colMeans(as.matrix(cluster_data[cluster_genes, ]))
    mean_expression_df <- data.frame(Sample = names(mean_expression), MeanExpression = mean_expression)
    plot_data <- merge(mean_expression_df, sample_mapping, by = "Sample", all.x = TRUE)
    mean_expression_data <- plot_data %>% group_by(Genotype) %>% summarize(MeanExpression = mean(MeanExpression))

    
    # Adjusting the factor levels in the 'Genotype' column based on the unique values in the data
    sample_mapping$Genotype <- factor(sample_mapping$Genotype, levels = unique(sample_mapping$Genotype))

# Then use this adjusted order in your plot generation
p_mean <- ggplot(plot_data, aes(x = Genotype, y = MeanExpression)) +
  geom_boxplot(aes(fill = Genotype), alpha = 0.5) +
  geom_smooth(data = mean_expression_data, aes(x = Genotype, y = MeanExpression, group = 1),
              method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  geom_point(data = mean_expression_data, aes(x = Genotype, y = MeanExpression),
             color = "black", size = 3) +
  labs(title = paste("Mean Expression Trend for Cluster", clust, "-", gene_count, "Genes"),
       x = "Genotype", y = "Mean Expression") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
  scale_x_discrete(limits = levels(sample_mapping$Genotype))  # Use the adjusted factor order

# Save the plot
tiff(paste0(opt$prefix, "_Cluster_", clust, "_", prefix_suffix, "_MeanExpression.tiff"), width = 10, height = 6, units = 'in', res = 300)
print(p_mean)
dev.off()

    
  }
}

######combined MeanExpression plot 
# Function to generate combined cluster Mean Expression trend plots
# Function to generate combined cluster Mean Expression trend plots
######combined MeanExpression plot 
# Function to generate combined cluster Mean Expression trend plots (CLUSTERS IN ORDER 1,2,3,...)
generate_combined_mean_expression_plots <- function(cluster_data, cluster_labels, prefix_suffix) {
  # Collect individual cluster plots
  cluster_plots <- list()
  
  # Loop through each cluster IN NUMERIC ORDER
  for (clust in sort(unique(cluster_labels))) {          # <-- CHANGE 1: sort() here
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    gene_count <- length(cluster_genes)
    
    # Skip empty clusters
    if (gene_count == 0) {
      next
    }

    # Ensure annotation follows the same order
    annotation_col <- data.frame(Genotype = factor(design_matrix$Genotype, levels = unique(design_matrix$Genotype)))
    rownames(annotation_col) <- colnames(cluster_data)

    # Generate Mean Expression data
    mean_expression <- colMeans(as.matrix(cluster_data[cluster_genes, ]))
    mean_expression_df <- data.frame(Sample = names(mean_expression), MeanExpression = mean_expression)
    plot_data <- merge(mean_expression_df, sample_mapping, by = "Sample", all.x = TRUE)
    mean_expression_data <- plot_data %>% group_by(Genotype) %>% summarize(MeanExpression = mean(MeanExpression))

    # Adjust the factor levels in the 'Genotype' column
    sample_mapping$Genotype <- factor(sample_mapping$Genotype, levels = unique(sample_mapping$Genotype))

    # Create the plot for each cluster
    p_mean <- ggplot(plot_data, aes(x = Genotype, y = MeanExpression)) +
      geom_boxplot(aes(fill = Genotype), alpha = 0.5) +
      geom_smooth(data = mean_expression_data, aes(x = Genotype, y = MeanExpression, group = 1),
                  method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
      geom_point(data = mean_expression_data, aes(x = Genotype, y = MeanExpression),
                 color = "black", size = 3) +
      labs(title = paste("Mean Expression Trend for Cluster", clust, "-", gene_count, "Genes"),
           x = "Genotype", y = "Mean Expression") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
      scale_x_discrete(limits = levels(sample_mapping$Genotype))
    
    # Store with character name (safe for list)
    cluster_plots[[as.character(clust)]] <- p_mean
  }

  # Check if we have any plots in the list
  if (length(cluster_plots) > 0) {
    # <-- CHANGE 2: Explicitly order the list by cluster number
    cluster_plots <- cluster_plots[order(as.numeric(names(cluster_plots)))]

    # Combine all cluster plots into one combined figure
    combined_plot <- plot_grid(plotlist = cluster_plots, ncol = 3, 
                               align = "hv", axis = "tblr")  # optional: better alignment
    
    # Save the combined plot
    tiff(paste0(opt$prefix, "_Combined_MeanExpression_", prefix_suffix, ".tiff"), 
         width = 15, height = 10, units = 'in', res = 300)
    print(combined_plot)
    dev.off()
    
    message("Combined plot saved with clusters ordered: ", 
            paste(names(cluster_plots), collapse = ", "))
  } else {
    warning("No valid clusters found for plotting.")
  }
}

# Call the function (unchanged)
generate_combined_mean_expression_plots(deg_data, filtered_clusters, "AfterFiltering")



# ==================================================================
# === CORRECTED HEATMAPS WITH AUTOMATIC Day0 → Day1 → Day3 → ... → Day70 ORDER === SORT genotype
# ==================================================================

# Function to sort Genotype levels chronologically (Day1 before Day10)
sort_genotype_levels <- function(genotypes) {
  day_nums <- as.numeric(sub("^Day", "", genotypes))
  genotypes[order(day_nums)]
}

# Apply correct order ONCE to design_matrix (affects ALL plots)
design_matrix$Genotype <- factor(design_matrix$Genotype,
                                 levels = sort_genotype_levels(unique(design_matrix$Genotype)))

# Rebuild annotation_col and sample_mapping with correct order
annotation_col <- data.frame(Genotype = design_matrix$Genotype)
rownames(annotation_col) <- design_matrix$Sample

sample_mapping <- data.frame(Sample = colnames(deg_data), 
                             Genotype = design_matrix$Genotype,
                             stringsAsFactors = FALSE)
sample_mapping$Genotype <- factor(sample_mapping$Genotype, 
                                  levels = levels(design_matrix$Genotype))


#####below here is before and after combined heatmap 
#####below here is before and after combined heatmap 
#####below here is before and after combined heatmap 
#####below here is before and after combined heatmap 
#####below here is before and after combined heatmap 
#####below here is before and after combined heatmap 

# Generate cluster heatmaps and trend plots before filtering
generate_cluster_plots(deg_data, mfuzz_result$cluster, "BeforeFiltering")

# Generate cluster heatmaps and trend plots after filtering
generate_cluster_plots(deg_data, filtered_clusters, "AfterFiltering")

# Generate heatmap for unfiltered clusters before and after filtering 
sorted_indices_before <- order(mfuzz_result$cluster)
scaled_data_before <- deg_data[sorted_indices_before, design_matrix$Sample]
annotation_row_before <- data.frame(Cluster = as.factor(mfuzz_result$cluster[sorted_indices_before]))
rownames(annotation_row_before) <- rownames(scaled_data_before)

gap_positions_before <- which(diff(as.numeric(annotation_row_before$Cluster)) != 0)

png(paste0(opt$prefix, "_AllClusters_BeforeFiltering.png"), width=2400, height=1800, res=300)
pheatmap(scaled_data_before, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,  # Preserve sample order
         show_rownames=FALSE, 
         annotation_row=annotation_row_before, 
         annotation_col=annotation_col, 
         gaps_row=gap_positions_before,  # Ensure proper gaps between clusters
         scale="row",  # Apply global row scaling
         color=colorRampPalette(c("blue", "white", "red"))(100),
         legend=FALSE)  # Blue to red scaling

dev.off()

# Reorder genes within each cluster while maintaining correct cluster annotations
sorted_indices <- order(filtered_clusters)
scaled_data <- deg_data[sorted_indices, design_matrix$Sample]
annotation_row <- data.frame(Cluster = as.factor(filtered_clusters[sorted_indices]))
rownames(annotation_row) <- rownames(scaled_data)
annotation_col <- data.frame(Genotype = design_matrix$Genotype)
rownames(annotation_col) <- colnames(scaled_data)

# Define gaps at boundaries between clusters
gap_positions <- which(diff(as.numeric(annotation_row$Cluster)) != 0)

# Generate heatmap with corrected cluster annotation, row scaling, and sample order maintained
png(paste0(opt$prefix, "_AllClusters_AfterFiltering4.png"), width=2400, height=1800, res=300)
pheatmap(scaled_data, 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,  # Preserve sample order
         show_rownames=FALSE, 
         annotation_row=annotation_row, 
         annotation_col=annotation_col, 
         gaps_row=gap_positions,  # Ensure proper gaps between clusters
         scale="row",  # Apply global row scaling
         annotation_names_col = FALSE,
         color=colorRampPalette(c("blue", "white", "red"))(100),
         main = "Combined Heatmap of all genes Clusters",
         legend=TRUE)  # Blue to red scaling

dev.off()



# ==================================================================
# === ALL CLUSTERS LINE TRAJECTORIES COMBINED (After Filtering) ===
# ==================================================================

generate_combined_line_trajectories <- function(cluster_data, cluster_labels, prefix_suffix) {
  # List to store individual trajectory plots
  trajectory_plots <- list()
  
  # Ensure Genotype is properly ordered (Day0 → Day1 → Day3 → ...)
  design_matrix$Genotype <- factor(design_matrix$Genotype,
                                   levels = sort_genotype_levels(unique(design_matrix$Genotype)))
  
  # Rebuild sample mapping with correct order
  sample_mapping <- data.frame(
    Sample = colnames(cluster_data),
    Genotype = design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)],
    NumericTime = as.numeric(gsub("^Day", "", design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)])),
    stringsAsFactors = FALSE
  )
  sample_mapping$Genotype <- factor(sample_mapping$Genotype, levels = levels(design_matrix$Genotype))
  
  # Create equidistant time index
  time_map <- sample_mapping %>%
    distinct(Genotype, NumericTime) %>%
    arrange(NumericTime) %>%
    mutate(TimeIndex = seq_len(n()))
  
  # Loop through clusters in numeric order
  for (clust in sort(unique(cluster_labels))) {
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    if (length(cluster_genes) == 0) next
    
    # Long format expression
    expr_long <- as.data.frame(t(cluster_data[cluster_genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
      left_join(sample_mapping, by = "Sample") %>%
      left_join(time_map, by = c("Genotype", "NumericTime"))
    
    # Average per gene per time point
    avg_per_gene <- expr_long %>%
      group_by(Gene, Genotype, TimeIndex) %>%
      summarise(Expression = mean(Expression), .groups = "drop")
    
    # Mean expression per gene (for color gradient)
    gene_means <- avg_per_gene %>%
      group_by(Gene) %>%
      summarise(GeneMean = mean(Expression), .groups = "drop")
    
    avg_per_gene <- avg_per_gene %>% left_join(gene_means, by = "Gene")
    
    # Plot
    p_traj <- ggplot(avg_per_gene, aes(x = TimeIndex, y = Expression, group = Gene, color = GeneMean)) +
      geom_line(alpha = 0.2, linewidth = 0.3) +
      scale_color_gradient(low = "#e0f3f8", high = "#0868ac", guide = "none") +
      scale_x_continuous(
        breaks = time_map$TimeIndex,
        labels = time_map$Genotype
      ) +
      labs(
        title = paste("Cluster", clust),
        subtitle = paste(length(cluster_genes), "genes"),
        x = NULL,
        y = "Standardized Expression (z-score)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    trajectory_plots[[as.character(clust)]] <- p_traj
  }
  
  # If no plots, exit
  if (length(trajectory_plots) == 0) {
    warning("No clusters with genes for line trajectories.")
    return()
  }
  
  # Order plots by cluster number
  trajectory_plots <- trajectory_plots[order(as.numeric(names(trajectory_plots)))]
  
  # Grid layout: square-like
  n <- length(trajectory_plots)
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  if (rows * cols < n) rows <- rows + 1
  
  # Combine
  combined <- cowplot::plot_grid(
    plotlist = trajectory_plots,
    nrow = rows, ncol = cols,
    align = "hv", axis = "tblr"
  )
  
  # Title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("All Clusters Line Trajectories – ", prefix_suffix, " (", n, " clusters)"),
      fontface = "bold", size = 16
    ) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  final_plot <- cowplot::plot_grid(title, combined, ncol = 1, rel_heights = c(0.06, 0.94))
  
  # Save PNG
  png(paste0(opt$prefix, "_AllClusters_LineTrajectories_Combined_", prefix_suffix, ".png"),
      width = max(12, cols * 4), height = max(8, rows * 3.5), units = "in", res = 300)
  print(final_plot)
  dev.off()
  
  # Save TIFF
  tiff(paste0(opt$prefix, "_AllClusters_LineTrajectories_Combined_", prefix_suffix, ".tiff"),
       width = max(12, cols * 4), height = max(8, rows * 3.5), units = "in", res = 300)
  print(final_plot)
  dev.off()
  
  message("Saved: ", opt$prefix, "_AllClusters_LineTrajectories_Combined_", prefix_suffix, ".png")
}

# ==================================================================
# === CALL THE FUNCTION (After Filtering) ==========================
# ==================================================================

generate_combined_line_trajectories(deg_data, filtered_clusters, "AfterFiltering")

# Optional: Before filtering
# generate_combined_line_trajectories(deg_data, mfuzz_result$cluster, "BeforeFiltering")


# ==================================================================
# === ALL CLUSTERS LINE TRAJECTORIES – COLORED BY MFUZZ MEMBERSHIP ===
# ==================================================================

generate_combined_line_trajectories <- function(cluster_data, cluster_labels, prefix_suffix) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  
  # === Automatically get membership from mfuzz_result (which already exists in your environment) ===
  if (!exists("mfuzz_result") || is.null(mfuzz_result$membership)) {
    stop("mfuzz_result$membership not found! Run Mfuzz clustering first.")
  }
  membership_matrix <- mfuzz_result$membership
  
  trajectory_plots <- list()
  
  # Ensure correct Day0 → Day1 → Day3 → ... order
  design_matrix$Genotype <- factor(design_matrix$Genotype,
                                   levels = sort_genotype_levels(unique(design_matrix$Genotype)))
  
  sample_mapping <- data.frame(
    Sample = colnames(cluster_data),
    Genotype = design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)],
    NumericTime = as.numeric(gsub("^Day", "", design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)])),
    stringsAsFactors = FALSE
  )
  sample_mapping$Genotype <- factor(sample_mapping$Genotype, levels = levels(design_matrix$Genotype))
  
  time_map <- sample_mapping %>%
    distinct(Genotype, NumericTime) %>%
    arrange(NumericTime) %>%
    mutate(TimeIndex = seq_len(n()))
  
  # Process each cluster
  for (clust in sort(unique(cluster_labels[cluster_labels != 0]))) {
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    if (length(cluster_genes) == 0) next
    
    # Get membership score for THIS cluster (column name = cluster number as character)
    gene_membership <- membership_matrix[cluster_genes, as.character(clust)]
    
    # Long format + join membership
    expr_long <- as.data.frame(t(cluster_data[cluster_genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
      left_join(sample_mapping, by = "Sample") %>%
      left_join(time_map, by = c("Genotype", "NumericTime")) %>%
      mutate(Membership = gene_membership[Gene])  # add membership
    
    # Average per gene per time point
    avg_per_gene <- expr_long %>%
      group_by(Gene, Genotype, TimeIndex, Membership) %>%
      summarise(Expression = mean(Expression), .groups = "drop")
    
    # Plot: color = membership score
    p <- ggplot(avg_per_gene, aes(x = TimeIndex, y = Expression, group = Gene, color = Membership)) +
      geom_line(alpha = 0.8, linewidth = 0.45) +
      scale_color_gradientn(
        colours = c("gray92", "#cce5ff", "#99ccff", "#3399ff", "#0066cc", "#003366"),
        values = scales::rescale(c(0, 0.3, 0.5, 0.7, 0.9, 1.0)),
        limits = c(0, 1),
        name = "Mfuzz\nMembership",
        breaks = c(0.2, 0.5, 0.8, 1.0)
      ) +
      scale_x_continuous(breaks = time_map$TimeIndex, labels = time_map$Genotype) +
      labs(
        title = paste("Cluster", clust),
        subtitle = paste(length(cluster_genes), "genes"),
        x = NULL,
        y = "Standardized Expression (z-score)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = if (clust == min(sort(unique(cluster_labels[cluster_labels != 0])))) "right" else "none"
      )
    
    trajectory_plots[[as.character(clust)]] <- p
  }
  
  if (length(trajectory_plots) == 0) {
    warning("No clusters to plot.")
    return(invisible(NULL))
  }
  
  # Order plots
  trajectory_plots <- trajectory_plots[order(as.numeric(names(trajectory_plots)))]
  
  # Extract legend from first plot
  legend <- cowplot::get_legend(trajectory_plots[[1]] + theme(legend.position = "right"))
  trajectory_plots <- lapply(trajectory_plots, function(p) p + theme(legend.position = "none"))
  
  # Arrange in grid
  n <- length(trajectory_plots)
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  
  grid <- cowplot::plot_grid(plotlist = trajectory_plots, nrow = rows, ncol = cols, align = "hv")
  with_legend <- cowplot::plot_grid(grid, legend, rel_widths = c(1, 0.2))
  
  # Title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Line Trajectories Colored by Mfuzz Membership – ", prefix_suffix, " (", n, " clusters)"),
      fontface = "bold", size = 18
    )
  
  final_plot <- cowplot::plot_grid(title, with_legend, ncol = 1, rel_heights = c(0.08, 0.92))
  
  # Save
  png(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".png"),
      width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  tiff(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".tiff"),
       width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  message("Saved membership-colored trajectories: ", 
          opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".png")
}

# ==================================================================
# === CALL THE FUNCTION (After Filtering) ==========================
# ==================================================================

generate_combined_line_trajectories(deg_data, filtered_clusters, "AfterFiltering")

# Optional: before filtering
generate_combined_line_trajectories(deg_data, mfuzz_result$cluster, "BeforeFiltering")


# # ==================================================================
# # === ALL CLUSTERS LINE TRAJECTORIES – COLORED BY MFUZZ MEMBERSHIP ===
# # ==================================================================

# generate_combined_line_trajectories <- function(cluster_data, cluster_labels, prefix_suffix) {
  
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(cowplot)
  
#   # === Automatically get membership from mfuzz_result (which already exists in your environment) ===
#   if (!exists("mfuzz_result") || is.null(mfuzz_result$membership)) {
#     stop("mfuzz_result$membership not found! Run Mfuzz clustering first.")
#   }
#   membership_matrix <- mfuzz_result$membership
  
#   trajectory_plots <- list()
  
#   # Ensure correct Day0 → Day1 → Day3 → ... order
#   # NOTE: sort_genotype_levels function must be defined elsewhere in your script
#   design_matrix$Genotype <- factor(design_matrix$Genotype,
#                                    levels = sort_genotype_levels(unique(design_matrix$Genotype)))
  
#   sample_mapping <- data.frame(
#     Sample = colnames(cluster_data),
#     Genotype = design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)],
#     NumericTime = as.numeric(gsub("^Day", "", design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)])),
#     stringsAsFactors = FALSE
#   )
#   sample_mapping$Genotype <- factor(sample_mapping$Genotype, levels = levels(design_matrix$Genotype))
  
#   time_map <- sample_mapping %>%
#     distinct(Genotype, NumericTime) %>%
#     arrange(NumericTime) %>%
#     mutate(TimeIndex = seq_len(n()))
  
#   # Process each cluster
#   for (clust in sort(unique(cluster_labels[cluster_labels != 0]))) {
#     cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
#     if (length(cluster_genes) == 0) next
    
#     # Get membership score for THIS cluster (column name = cluster number as character)
#     gene_membership <- membership_matrix[cluster_genes, as.character(clust)]
    
#     # Long format + join membership
#     expr_long <- as.data.frame(t(cluster_data[cluster_genes, , drop = FALSE])) %>%
#       tibble::rownames_to_column("Sample") %>%
#       pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
#       left_join(sample_mapping, by = "Sample") %>%
#       left_join(time_map, by = c("Genotype", "NumericTime")) %>%
#       mutate(Membership = gene_membership[Gene])  # add membership
    
#     # Average per gene per time point
#     avg_per_gene <- expr_long %>%
#       group_by(Gene, Genotype, TimeIndex, Membership) %>%
#       summarise(Expression = mean(Expression), .groups = "drop")
    
#     # Plot: color = membership score
#     p <- ggplot(avg_per_gene, aes(x = TimeIndex, y = Expression, group = Gene, color = Membership)) +
#       geom_line(alpha = 0.8, linewidth = 0.45) +
#       scale_color_gradientn(
#         colours = c("gray92", "#cce5ff", "#99ccff", "#3399ff", "#0066cc", "#003366"),
#         values = scales::rescale(c(0, 0.3, 0.5, 0.7, 0.9, 1.0)),
#         limits = c(0, 1),
#         name = "Mfuzz\nMembership",
#         breaks = c(0.2, 0.5, 0.8, 1.0)
#       ) +
#       scale_x_continuous(breaks = time_map$TimeIndex, labels = time_map$Genotype) +
#       labs(
#         title = paste("Cluster", clust),
#         subtitle = paste(length(cluster_genes), "genes"),
#         x = NULL,
#         y = "Standardized Expression (z-score)"
#       ) +
#       theme_classic() +
#       theme(
#         plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#         plot.subtitle = element_text(size = 11, hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#         legend.position = if (clust == min(sort(unique(cluster_labels[cluster_labels != 0])))) "right" else "none"
#       )
    
#     trajectory_plots[[as.character(clust)]] <- p
#   }
  
#   if (length(trajectory_plots) == 0) {
#     warning("No clusters to plot.")
#     return(invisible(NULL))
#   }
  
#   # Order plots
#   trajectory_plots <- trajectory_plots[order(as.numeric(names(trajectory_plots)))]
  
#   # Extract legend from first plot
#   legend <- cowplot::get_legend(trajectory_plots[[1]] + theme(legend.position = "right"))
#   trajectory_plots <- lapply(trajectory_plots, function(p) p + theme(legend.position = "none"))
  
#   # Arrange in grid
#   n <- length(trajectory_plots)
#   cols <- ceiling(sqrt(n))
#   rows <- ceiling(n / cols)
  
#   grid <- cowplot::plot_grid(plotlist = trajectory_plots, nrow = rows, ncol = cols, align = "hv")
#   with_legend <- cowplot::plot_grid(grid, legend, rel_widths = c(1, 0.2))
  
#   # Title
#   title <- cowplot::ggdraw() +
#     cowplot::draw_label(
#       paste0("Line Trajectories Colored by Mfuzz Membership – ", prefix_suffix, " (", n, " clusters)"),
#       fontface = "bold", size = 18
#     )
  
#   final_plot <- cowplot::plot_grid(title, with_legend, ncol = 1, rel_heights = c(0.08, 0.92))
  
#   # Save
#   # NOTE: opt$prefix variable must be defined elsewhere
#   png(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorPlor_", prefix_suffix, ".png"),
#       width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
#   print(final_plot)
#   dev.off()
  
#   tiff(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorPlot_", prefix_suffix, ".tiff"),
#        width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
#   print(final_plot)
#   dev.off()
  
#   message("Saved membership-colored trajectories: ", 
#           opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".png")
# }

# # ==================================================================
# # === CALL THE FUNCTION (After Filtering) ==========================
# # ==================================================================

# generate_combined_line_trajectories(deg_data, filtered_clusters, "AfterFiltering")

# # Optional: before filtering
# generate_combined_line_trajectories(deg_data, mfuzz_result$cluster, "BeforeFiltering")


###perflexity 

# ==================================================================
# === ALL CLUSTERS LINE TRAJECTORIES – COLORED BY MFUZZ MEMBERSHIP ===
# ==================================================================

generate_combined_line_trajectories <- function(cluster_data, cluster_labels, prefix_suffix) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  
  # === Automatically get membership from mfuzz_result ===
  if (!exists("mfuzz_result") || is.null(mfuzz_result$membership)) {
    stop("mfuzz_result$membership not found! Run Mfuzz clustering first.")
  }
  membership_matrix <- mfuzz_result$membership
  
  trajectory_plots <- list()
  
  # Ensure correct Day0 → Day1 → Day3 → ... order
  design_matrix$Genotype <- factor(
    design_matrix$Genotype,
    levels = sort_genotype_levels(unique(design_matrix$Genotype))
  )
  
  sample_mapping <- data.frame(
    Sample      = colnames(cluster_data),
    Genotype    = design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)],
    NumericTime = as.numeric(gsub("^Day", "", 
                                  design_matrix$Genotype[match(colnames(cluster_data), 
                                                               design_matrix$Sample)])),
    stringsAsFactors = FALSE
  )
  sample_mapping$Genotype <- factor(sample_mapping$Genotype,
                                    levels = levels(design_matrix$Genotype))
  
  time_map <- sample_mapping %>%
    distinct(Genotype, NumericTime) %>%
    arrange(NumericTime) %>%
    mutate(TimeIndex = seq_len(n()))
  
  # Process each cluster
  for (clust in sort(unique(cluster_labels[cluster_labels != 0]))) {
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    if (length(cluster_genes) == 0) next
    
    # Membership score for this cluster
    gene_membership <- membership_matrix[cluster_genes, as.character(clust)]
    
    # Long format + join membership
    expr_long <- as.data.frame(t(cluster_data[cluster_genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
      left_join(sample_mapping, by = "Sample") %>%
      left_join(time_map, by = c("Genotype", "NumericTime")) %>%
      mutate(Membership = gene_membership[Gene])
    
    # Average per gene per time point
    avg_per_gene <- expr_long %>%
      group_by(Gene, Genotype, TimeIndex, Membership) %>%
      summarise(Expression = mean(Expression), .groups = "drop")
    
    ######Plot: color = membership score, Spectral‑like palette
    p <- ggplot(avg_per_gene,
                aes(x = TimeIndex, y = Expression, group = Gene, color = Membership)) +
      geom_line(alpha = 0.8, linewidth = 0.45) +
      scale_color_gradientn(
        colours = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
                    "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
                    "#F46D43", "#D73027", "#A50026"),
        values  = scales::rescale(seq(0, 1, length.out = 11)),
        limits  = c(0, 1),

      name    = "Mfuzz\nMembership"
      ) +
      scale_x_continuous(
        breaks = time_map$TimeIndex,
        labels = time_map$Genotype
      ) +
      labs(
        title    = paste("Cluster", clust),
        subtitle = paste(length(cluster_genes), "genes"),
        x = NULL,
        y = "Standardized Expression (z-score)"
      ) +
      theme_classic() +
      theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = if (clust == min(sort(unique(cluster_labels[cluster_labels != 0]))))
          "right" else "none"
      )
    
    trajectory_plots[[as.character(clust)]] <- p
  }
  
  if (length(trajectory_plots) == 0) {
    warning("No clusters to plot.")
    return(invisible(NULL))
  }
  
  # Order plots
  trajectory_plots <- trajectory_plots[order(as.numeric(names(trajectory_plots)))]
  
  # Extract legend from first plot
  legend <- cowplot::get_legend(trajectory_plots[[1]] + theme(legend.position = "right"))
  trajectory_plots <- lapply(trajectory_plots, function(p) p + theme(legend.position = "none"))
  
  # Arrange in grid (same panel structure as before)
  n    <- length(trajectory_plots)
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  
  grid <- cowplot::plot_grid(
    plotlist = trajectory_plots,
    nrow = rows, ncol = cols,
    align = "hv"
  )
  with_legend <- cowplot::plot_grid(grid, legend, rel_widths = c(1, 0.2))
  
  # Title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Line Trajectories Colored by Mfuzz Membership – ",
             prefix_suffix, " (", n, " clusters)"),
      fontface = "bold", size = 18
    )
  
  final_plot <- cowplot::plot_grid(
    title, with_legend,
    ncol = 1, rel_heights = c(0.08, 0.92)
  )
  
  # Save
  png(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorFacetRed_", prefix_suffix, ".png"),
      width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  tiff(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorFacetRed_", prefix_suffix, ".tiff"),
       width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  message("Saved membership-colored trajectories: ",
          opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".png")
}

# ==================================================================
# === CALL THE FUNCTION (After Filtering) ==========================
# ==================================================================

generate_combined_line_trajectories(deg_data, filtered_clusters, "AfterFiltering")

# Optional: before filtering
generate_combined_line_trajectories(deg_data, mfuzz_result$cluster, "BeforeFiltering")




###rainbow colors 

# ==================================================================
# === ALL CLUSTERS LINE TRAJECTORIES – COLORED BY MFUZZ MEMBERSHIP ===
# ==================================================================

generate_combined_line_trajectories <- function(cluster_data, cluster_labels, prefix_suffix) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  
  # === Automatically get membership from mfuzz_result ===
  if (!exists("mfuzz_result") || is.null(mfuzz_result$membership)) {
    stop("mfuzz_result$membership not found! Run Mfuzz clustering first.")
  }
  membership_matrix <- mfuzz_result$membership
  
  trajectory_plots <- list()
  
  # Ensure correct Day0 → Day1 → Day3 → ... order
  design_matrix$Genotype <- factor(
    design_matrix$Genotype,
    levels = sort_genotype_levels(unique(design_matrix$Genotype))
  )
  
  sample_mapping <- data.frame(
    Sample      = colnames(cluster_data),
    Genotype    = design_matrix$Genotype[match(colnames(cluster_data), design_matrix$Sample)],
    NumericTime = as.numeric(gsub("^Day", "", 
                                  design_matrix$Genotype[match(colnames(cluster_data), 
                                                               design_matrix$Sample)])),
    stringsAsFactors = FALSE
  )
  sample_mapping$Genotype <- factor(sample_mapping$Genotype,
                                    levels = levels(design_matrix$Genotype))
  
  time_map <- sample_mapping %>%
    distinct(Genotype, NumericTime) %>%
    arrange(NumericTime) %>%
    mutate(TimeIndex = seq_len(n()))
  
  # Process each cluster
  for (clust in sort(unique(cluster_labels[cluster_labels != 0]))) {
    cluster_genes <- rownames(cluster_data)[cluster_labels == clust]
    if (length(cluster_genes) == 0) next
    
    # Membership score for this cluster
    gene_membership <- membership_matrix[cluster_genes, as.character(clust)]
    
    # Long format + join membership
    expr_long <- as.data.frame(t(cluster_data[cluster_genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
      left_join(sample_mapping, by = "Sample") %>%
      left_join(time_map, by = c("Genotype", "NumericTime")) %>%
      mutate(Membership = gene_membership[Gene])
    
    # Average per gene per time point
    avg_per_gene <- expr_long %>%
      group_by(Gene, Genotype, TimeIndex, Membership) %>%
      summarise(Expression = mean(Expression), .groups = "drop")
    

    ##different colors Expression (rainbow color) 
      p <- ggplot(avg_per_gene,
            aes(x = TimeIndex,
                y = Expression,
                group = Gene,
                color = Expression)) +    # color by expression
      geom_line(aes(alpha = Membership),     # optional: encode membership here
            linewidth = 0.45) +
      scale_color_gradientn(
      colours = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
                "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
                "#F46D43", "#D73027", "#A50026"),
      values  = scales::rescale(seq(0, 1, length.out = 11)),
      limits  = range(avg_per_gene$Expression, na.rm = TRUE),

      name    = "Mfuzz\nExpression"
      ) +
      scale_x_continuous(
        breaks = time_map$TimeIndex,
        labels = time_map$Genotype
      ) +
      labs(
        title    = paste("Cluster", clust),
        subtitle = paste(length(cluster_genes), "genes"),
        x = NULL,
        y = "Standardized Expression (z-score)"
      ) +
      theme_classic() +
      theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = if (clust == min(sort(unique(cluster_labels[cluster_labels != 0]))))
          "right" else "none"
      )
    
    trajectory_plots[[as.character(clust)]] <- p
  }
  
  if (length(trajectory_plots) == 0) {
    warning("No clusters to plot.")
    return(invisible(NULL))
  }
  
  # Order plots
  trajectory_plots <- trajectory_plots[order(as.numeric(names(trajectory_plots)))]
  
  # Extract legend from first plot
  legend <- cowplot::get_legend(trajectory_plots[[1]] + theme(legend.position = "right"))
  trajectory_plots <- lapply(trajectory_plots, function(p) p + theme(legend.position = "none"))
  
  # Arrange in grid (same panel structure as before)
  n    <- length(trajectory_plots)
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  
  grid <- cowplot::plot_grid(
    plotlist = trajectory_plots,
    nrow = rows, ncol = cols,
    align = "hv"
  )
  with_legend <- cowplot::plot_grid(grid, legend, rel_widths = c(1, 0.2))
  
  # Title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("Line Trajectories Colored by Mfuzz Membership – ",
             prefix_suffix, " (", n, " clusters)"),
      fontface = "bold", size = 18
    )
  
  final_plot <- cowplot::plot_grid(
    title, with_legend,
    ncol = 1, rel_heights = c(0.08, 0.92)
  )
  
  # Save
  png(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorFacet_", prefix_suffix, ".png"),
      width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  tiff(paste0(opt$prefix, "_AllClusters_LineTrajectories_MembershipColorFacet_", prefix_suffix, ".tiff"),
       width = max(14, cols * 4.5), height = max(9, rows * 4), units = "in", res = 320)
  print(final_plot)
  dev.off()
  
  message("Saved membership-colored trajectories: ",
          opt$prefix, "_AllClusters_LineTrajectories_MembershipColor_", prefix_suffix, ".png")
}

# ==================================================================
# === CALL THE FUNCTION (After Filtering) ==========================
# ==================================================================

generate_combined_line_trajectories(deg_data, filtered_clusters, "AfterFiltering")

# Optional: before filtering
generate_combined_line_trajectories(deg_data, mfuzz_result$cluster, "BeforeFiltering")

# Perform PCA
pca_res <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

# Create PCA dataframe
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Cluster = as.factor(filtered_clusters))

# Plot PCA
png(paste0(opt$prefix, "_PCA.png"), width=2400, height=1800, res=300)

ggplot(na.omit(pca_df), aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA - Gene Clusters")
dev.off()


# Perform UMAP
set.seed(123)
umap_res <- umap(scaled_data)
umap_df <- data.frame(UMAP1 = umap_res$layout[,1], UMAP2 = umap_res$layout[,2], Cluster = filtered_clusters)

# Plot UMAP
png(paste0(opt$prefix, "_UMAP.png"), width=2400, height=1800, res=300)
ggplot(na.omit(umap_df), aes(x=UMAP1, y=UMAP2, color=as.factor(Cluster))) + geom_point(alpha=0.7) + theme_minimal() + labs(title="UMAP - Gene Clusters")
dev.off()

# Perform t-SNE (perplexity is dataset dependent, adjust if necessary)

# Determine the number of samples
num_samples <- ncol(deg_data)

set.seed(123)
# Set perplexity dynamically (max value is num_samples/3)
perplexity_value <- min(30, floor(num_samples / 3))

# Perform t-SNE
set.seed(123)
tsne_res <- Rtsne(scale(deg_data), perplexity = perplexity_value, theta = 0.5, dims = 2)

# Create t-SNE dataframe
tsne_df <- data.frame(tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2], Cluster = as.factor(filtered_clusters))

# Plot t-SNE
png(paste0(opt$prefix, "_tSNE.png"), width=2400, height=1800, res=300)
ggplot(na.omit(tsne_df), aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE - Gene Clusters")
dev.off()
