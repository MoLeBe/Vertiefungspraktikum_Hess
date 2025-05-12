# ============================================
# insert_genome_seg
# Input:
#   - inp: Summarized_Experiment object containing RNA-seq data with genomic coordinates in rowRanges
#   - gff3_file_path: path to GFF3 annotation file
# Output: Modified SummarizedExperiment object with added 
#         gene-level segment annotations (gene_segment_ID and gene_name) in rowRanges
# ============================================

insert_genome_seg <- function(inp, gff3_file_path) {
  # Load and preprocess GFF3 data
  gff3_data <- gff3_preprocess(gff3_file_path) %>%
    arrange(start) %>%
    mutate(gene_segment_ID = paste0("G_", row_number()))  # Assign unique IDs
  
  # Convert to GRanges object
  gff3_gr <- GRanges(
    seqnames = "chr",  # Adjust if chromosome info is available
    ranges = IRanges(start = gff3_data$start, end = gff3_data$end),
    strand = gff3_data$strand,
    gene_locus = gff3_data$gene_locus
  )
  
  # Ensure unique gene_segment_IDs by gene_locus and strand
  gene_ids <- gff3_data %>%
    distinct(strand, gene_locus) %>%
    group_by(strand) %>%
    mutate(gene_segment_ID = paste0("G_", row_number())) %>%
    ungroup()
  
  # Annotate GRanges with IDs
  gff3_gr$gene_segment_ID <- gene_ids$gene_segment_ID[match(gff3_gr$gene_locus, gene_ids$gene_locus)]
  
  # Annotate rowRanges of input object
  rr <- rowRanges(inp)
  hits <- findOverlaps(rr, gff3_gr, ignore.strand = FALSE)
  rr$gene_segment_ID <- NA_character_
  rr$gene_segment_ID[queryHits(hits)] <- gff3_gr$gene_segment_ID[subjectHits(hits)]
  rr$gene_name <- NA_character_
  rr$gene_name[queryHits(hits)] <- gff3_gr$gene_locus[subjectHits(hits)]
  
  # Update and return
  rowRanges(inp) <- rr
  return(inp)
}


preprocess_comp_se_gene <- function(comp_se, strand_type) {
  rowRanges(comp_se) %>%
    as.data.frame() %>%
    filter(if (strand_type == "both") TRUE else strand == strand_type) %>%
    select(
      half_life_dif, log2_half_life_fold_change, log2_intensity_fold_change,
      position, decay_constant_1, decay_constant_2,
      synthesis_rate_1, synthesis_rate_2, ID, intensity, intensity_2,
      strand, gene_segment_ID, gene_name
    )
}

# ============================================
# segment_gene_comp
# Input:
#   - se_cond_1: SummarizedExperiment for condition 1
#   - se_cond_2: SummarizedExperiment for condition 
#   - gff3_file_path: path to GFF3 annotation file
# Output: SummarizedExperiment with added difference columns (e.g., half-life difference, fold-changes).
# ============================================


segment_gene_comp <- function(se_cond_1, se_cond_2) {
  se_cond_1 <- decay_c_synth(se_cond_1)
  se_cond_2 <- decay_c_synth(se_cond_2)
  
  rr1 <- rowRanges(se_cond_1)
  rr2 <- rowRanges(se_cond_2)
  
  # Compute differential metrics
  rr1$half_life_dif <- rr2$half_life - rr1$half_life
  rr1$log2_half_life_fold_change <- log2(rr2$half_life / rr1$half_life)
  rr1$decay_constant_1 <- rr1$decay_constant
  rr1$decay_constant_2 <- rr2$decay_constant
  rr1$synthesis_rate_1 <- rr1$synthesis_rate
  rr1$synthesis_rate_2 <- rr2$synthesis_rate
  rr1$intensity_2 <- rr2$intensity
  rr1$log2_intensity_fold_change <- log2(rr1$intensity_2 / rr1$intensity)
  
  rowRanges(se_cond_1) <- rr1
  se_cond_1 <- se_cond_1[!is.na(rr1$gene_segment_ID), ]
  return(se_cond_1)
}



group_by_seg_gene <- function(df) {
  # Remove rows where gene_segment_ID is missing
  df <- df %>%
    filter(!is.na(gene_segment_ID)) %>%
    group_by(gene_segment_ID) %>%
    mutate(
      across(
        c(half_life_dif, log2_half_life_fold_change,
          log2_intensity_fold_change, 
          decay_constant_1, 
          decay_constant_2, synthesis_rate_1, 
          synthesis_rate_2, intensity, 
          intensity_2),
        ~ detect_outliers_cleaned(.x)$cleaned_data,
        .names = "{.col}_clean"
      )
    ) %>%
    summarise(
      across(ends_with("_clean"), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
      position_1 = min(position, na.rm = TRUE),
      position_2 = max(position, na.rm = TRUE),
      mean_log_decay_constant_fold_change = log(mean(decay_constant_1_clean, na.rm = TRUE) / mean(decay_constant_2_clean, na.rm = TRUE)),
      mean_log_synthesis_rate_fold_change = log(mean(synthesis_rate_2_clean, na.rm = TRUE) / mean(synthesis_rate_1_clean, na.rm = TRUE)),
      mean_intensity = mean(intensity, na.rm = TRUE),
      mean_intensity_2 = mean(intensity_2, na.rm = TRUE),
      IDs = paste(unique(ID), collapse = ", "),
      gene_name = first(gene_name)
    ) %>%
    ungroup()
  
  return(df)
}
get_dif_se_gene <- function(se_cond_1, se_cond_2, gff3_file_path) {
  start_time <- Sys.time()
  
  se_cond_1_gene <- insert_genome_seg(se_cond_1, gff3_file_path)
  comp_se <- segment_gene_comp(se_cond_1_gene, se_cond_2)
  full_tmp_df <- preprocess_comp_se_gene(comp_se, "both")
  grouped_full_comp_df <- group_by_seg_gene(full_tmp_df)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(grouped_full_comp_df)
}



round_up <- function(x, precision = 0.1) {
  ceiling(x / precision) * precision
}
# ============================================
# hl_int_plot_gene
# Input:
# inp: A data frame containing the following columns:
#   - log2_intensity_fold_change_clean_mean: The mean log2 fold change in intensity
#   - log2_half_life_fold_change_clean_mean: The mean log2 fold change in half-life
#   - gene_name: The gene names
#   - half_life_dif_clean_mean: The mean difference in half-life
#   - log2_intensity_fold_change_clean_mean: The clean mean log2 fold change of intensity
# Output:
#   A PDF file with a plot showing the relationship between log2 intensity and log2 half-life fold changes
#   for various genes, with classification based on the ratio of these values.
# ============================================

hl_int_plot_gene <- function(inp, name) {
  delta = 0.5
  custom_colors <- c(
    "high d half_life low d_int" = "indianred",     # Red
    "low d half_life high d_int" = "steelblue",     # Blue
    "balanced" = "chartreuse4",                     # Green
    "overcompensation" = "gold",
    "other" = "darkgrey"                            # Grey
  )
  legend_labels <- c(
    "high d half_life low d_int" = expression("high" ~ Delta * "Half-Life, low" ~ Delta * "Intensity"),
    "low d half_life high d_int" = expression("low" ~ Delta * "Half-Life, high" ~ Delta * "Intensity"),
    "balanced" = "Balanced",
    "overcompensation" = "Overcompensation",
    "other" = "Other"
  )
  plot_data <- inp %>%
    mutate(
      label = if_else(
        (log2_intensity_fold_change_clean_mean >= 0.5 | log2_intensity_fold_change_clean_mean <= -0.5) & 
          (log2_half_life_fold_change_clean_mean >= 0.5 | log2_half_life_fold_change_clean_mean <= -0.5), 
        TRUE, 
        FALSE
      ),
      gene_name = str_replace_all(gene_name, ", NA|NA,", ""),
      ratio_value = half_life_dif_clean_mean / log2_intensity_fold_change_clean_mean,
      ratio = case_when(
        # Classify as "overcompensation" when close to y = -x line
        abs(log2_half_life_fold_change_clean_mean + log2_intensity_fold_change_clean_mean) < delta ~ "overcompensation",
        
        # Classify as "balanced" when close to y = x line
        abs(log2_half_life_fold_change_clean_mean - log2_intensity_fold_change_clean_mean) < delta ~ "balanced",
        
        # High d half-life, low d intensity
        abs(ratio_value) >= 3 ~ "high d half_life low d_int",
        
        # Low d half-life, high d intensity
        abs(ratio_value) <= 0.15 ~ "low d half_life high d_int",
        
        # Default to "other"
        TRUE ~ "other"
      )
    )
  # Extract top genes for high d_half_life low d_int
  top_high_hl_positive <- plot_data %>%
    filter(ratio == "high d half_life low d_int", log2_half_life_fold_change_clean_mean > 0) %>%
    slice_max(order_by = log2_half_life_fold_change_clean_mean, n = 5)
  
  top_high_hl_negative <- plot_data %>%
    filter(ratio == "high d half_life low d_int", log2_half_life_fold_change_clean_mean < 0) %>%
    slice_min(order_by = log2_half_life_fold_change_clean_mean, n = 5)
  
  # Extract top genes for low d_half_life high d_int
  top_low_hl_positive <- plot_data %>%
    filter(ratio == "low d half_life high d_int", log2_intensity_fold_change_clean_mean > 0) %>%
    slice_max(order_by = log2_intensity_fold_change_clean_mean, n = 5)
  
  top_low_hl_negative <- plot_data %>%
    filter(ratio == "low d half_life high d_int", log2_intensity_fold_change_clean_mean < 0) %>%
    slice_min(order_by = log2_intensity_fold_change_clean_mean, n = 5)
  
  # Combine all top genes
  top_high_hl <- bind_rows(top_high_hl_positive, top_high_hl_negative)
  top_low_hl <- bind_rows(top_low_hl_positive, top_low_hl_negative)
  top_genes <- bind_rows(top_high_hl, top_low_hl)
  
  
  # Determine plot dimensions
  overall_max <- round_up(max(
    pmax(
      abs(plot_data$log2_half_life_fold_change_clean_mean),
      abs(plot_data$log2_intensity_fold_change_clean_mean),
      na.rm = TRUE
    ),
    na.rm = TRUE
  ))
  x_dim <- 10
  y_dim <- 10
  
  # Create the base plot
  pdf(paste0("~/Desktop/R_plots/", name, ".pdf"), width = 10, height = 10)
  
  full_plot <- ggplot(data = plot_data, 
                      aes(x = log2_intensity_fold_change_clean_mean, 
                          y = log2_half_life_fold_change_clean_mean)) +
    geom_point(data = filter(plot_data, ratio != "other"), aes(color = ratio), size = 1.5) +
    geom_point(data = filter(plot_data, ratio == "other"), aes(color = ratio), size = 1) +
    geom_text(data = filter(plot_data, ratio == "balanced" & abs(log2_intensity_fold_change_clean_mean) > 0.8 & abs(half_life_dif_clean_mean) > 0.8),
              aes(label = gene_name),
              check_overlap = TRUE,
              nudge_y = 0.25,
              size = 3) +
    geom_text(data = top_genes,
              aes(label = gene_name),
              color = "black",
              check_overlap = TRUE,
              nudge_y = 0.3,
              size = 3) +  # Highlight top genes
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_abline(intercept = 0, slope = -1,linetype = "dashed", color = "black", alpha = 0.3 ) +
    geom_abline(intercept = 0, slope = 1,linetype = "dashed", color = "black" , alpha = 0.3) +
    xlim(-x_dim, x_dim) +
    ylim(-y_dim, y_dim) +
    xlab("Mean log2 intensity Fold Change") + 
    ylab("Mean log2 half-life Fold Change") +
    theme_bw() +
    scale_color_manual(values = custom_colors, labels = legend_labels) +
    labs(color = "Ratio") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12),      # Adjust size for legend labels
          legend.title = element_text(size = 14, face = "bold"),
          axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"),
          plot.title = element_text(size=14,face="bold")) +
    ggtitle("Differential Expression")
  
  # Print the plot
  print(full_plot)
  
  # Close the PDF device
  dev.off()
}


