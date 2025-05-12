
# ============================================
# hl_int_plot
# Input:
#   - inp: A dataframe containing differential expression data with columns:
#     - `log2_intensity_fold_change_clean_mean`: Mean log2 fold change of intensity
#     - `log2_half_life_fold_change_clean_mean`: Mean log2 fold change of half-life
#   - name: A string for naming the output plot file (PDF)
# Output:
#   - A ggplot object representing the plot of differential expression
#   - The plot is saved as a PDF file at "~/Desktop/R_plots/{name}.pdf"
# ============================================
hl_int_plot <- function(inp, name) {
  custom_colors <- c(
    "high d half_life low d_int" = "indianred",     # Red
    "low d half_life high d_int" = "steelblue",     # Blue
    "balanced" = "chartreuse4",                     # Green
    "overcompensation" = "gold",                    # Gold
    "other" = "darkgrey",                           # Grey
    "extremely high d_int" = "turquoise"            # Turquoise for extremely high delta intensity fold change
  )
  
  legend_labels <- c(
    "high d half_life low d_int" = expression("high" ~ Delta * "Half-Life, low" ~ Delta * "Intensity"),
    "low d half_life high d_int" = expression("low" ~ Delta * "Half-Life, high" ~ Delta * "Intensity"),
    "balanced" = "Balanced",
    "overcompensation" = "Overcompensation",
    "other" = "Other",
    "extremely high d_int" = expression("extremely high" ~ Delta * "Intensity")  # Label for the new category
  )
  
  delta = 0.8
  delta2 = 0.2
  
  plot_data <- inp$full_differential_expression %>%
    mutate(
      label = if_else(
        (log2_intensity_fold_change_clean_mean >= 0.5 | log2_intensity_fold_change_clean_mean <= -0.5) & 
          (log2_half_life_fold_change_clean_mean >= 0.5 | log2_half_life_fold_change_clean_mean <= -0.5), 
        TRUE, 
        FALSE
      ),
      gene = str_replace_all(gene, ", NA|NA,", ""),
      ratio_value = log2_half_life_fold_change_clean_mean / log2_intensity_fold_change_clean_mean,
      ratio = case_when(
        # Classify as "overcompensation" when close to y = -x line
        abs(log2_half_life_fold_change_clean_mean + log2_intensity_fold_change_clean_mean) < delta & 
          abs(log2_half_life_fold_change_clean_mean) > 0.4 & 
          abs(log2_intensity_fold_change_clean_mean) > 0.4 ~ "overcompensation",
        
        # Classify as "balanced" when close to y = x line
        abs(log2_half_life_fold_change_clean_mean - log2_intensity_fold_change_clean_mean) < delta & 
          abs(log2_half_life_fold_change_clean_mean) > 0.4 & 
          abs(log2_intensity_fold_change_clean_mean) > 0.4 ~ "balanced",
        
        # Classify as "high d half_life low d_int" when close to the y-axis (high half-life, low intensity)
        abs(log2_intensity_fold_change_clean_mean) < delta2 & 
          abs(log2_half_life_fold_change_clean_mean) > 0.4 ~ "high d half_life low d_int",
        
        # Classify as "low d half_life high d_int" when close to the x-axis (low half-life, high intensity)
        abs(log2_half_life_fold_change_clean_mean) < delta2 & 
          abs(log2_intensity_fold_change_clean_mean) > 0.4 ~ "low d half_life high d_int",
        
        # Classify as "extremely high d_int" for extremely high intensity fold change
        log2_intensity_fold_change_clean_mean > 3.5 ~ "extremely high d_int",
        
        # Default to "other"
        TRUE ~ "other"
      )
    )
  
  # Extract top genes for high d_half_life low d_int
  top_high_hl_positive <- plot_data %>%
    filter(
      ratio == "high d half_life low d_int",
      log2_half_life_fold_change_clean_mean > 0,
      !is.na(gene),
      gene != "NA"
    ) %>%
    slice_max(order_by = log2_half_life_fold_change_clean_mean, n = 10)
  
  top_high_hl_negative <- plot_data %>%
    filter(
      ratio == "high d half_life low d_int",
      log2_half_life_fold_change_clean_mean < 0,
      !is.na(gene),
      gene != "NA"
    ) %>%
    slice_min(order_by = log2_half_life_fold_change_clean_mean, n = 10)
  
  top_low_hl_positive <- plot_data %>%
    filter(
      ratio == "low d half_life high d_int",
      log2_intensity_fold_change_clean_mean > 0,
      !is.na(gene),
      gene != "NA"
    ) %>%
    slice_max(order_by = log2_intensity_fold_change_clean_mean, n = 10)
  
  top_low_hl_negative <- plot_data %>%
    filter(
      ratio == "low d half_life high d_int",
      log2_intensity_fold_change_clean_mean < 0,
      !is.na(gene),
      gene != "NA"
    ) %>%
    slice_min(order_by = log2_intensity_fold_change_clean_mean, n = 10)
  
  
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
  x_dim <- 7.5
  y_dim <- 7
  
  # Create the base plot
  pdf(paste0("~/Desktop/R_plots/", name, ".pdf"), width = 10, height = 10)
  full_plot <- ggplot(data = plot_data, 
                      aes(x = log2_intensity_fold_change_clean_mean, 
                          y = log2_half_life_fold_change_clean_mean)) +
    geom_point(data = filter(plot_data, ratio != "other"), aes(color = ratio), size = 1.5) +
    geom_point(data = filter(plot_data, ratio == "other"), aes(color = ratio), size = 1) +
    geom_text(data = filter(plot_data, 
                            ratio %in% c("balanced", "overcompensation", "extremely high d_int")),
              aes(label = gene),
              check_overlap = TRUE,
              nudge_y = 0.25,
              size = 3) +
    geom_text(data = top_genes,
              aes(label = gene),
              color = "black",
              check_overlap = TRUE,
              nudge_y = 0.3,
              size = 3) +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black", alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.3) +
    xlim(-4, 7.5) +
    ylim(-y_dim, y_dim) +
    xlab("Mean log2 intensity Fold Change") + 
    ylab("Mean log2 half life Fold Change") +
    theme_bw() +
    scale_color_manual(values = custom_colors, labels = legend_labels,
                       breaks = c("high d half_life low d_int", "low d half_life high d_int", "balanced", "overcompensation", "other", "extremely high d_int")) +
    labs(color = "Relationship") +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),      
          legend.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 14, face = "bold")) +
    ggtitle("Differential Expression")
  
  # Print the plot
  print(full_plot)
  
  # Close the PDF device
  dev.off()
  
  return(full_plot)
}

