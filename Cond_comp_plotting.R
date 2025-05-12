
# Add Roboto font (download if needed)
font_add_google("Roboto", "Roboto")

# Enable showtext
showtext_auto()

# Outlier definition function
get_hlines <- function(df, y_col) {
  y_col <- as.character(substitute(y_col))
  column_data <- df[[y_col]]
  
  # --- Outlier detection (same logic as in process_segment) ---
  mad_value <- mad(column_data, constant = 1.4826, na.rm = TRUE)
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- 1
  }
  median_value <- median(column_data, na.rm = TRUE)
  
  modified_z_score <- 0.6745 * (column_data - median_value) / mad_value
  non_outlier_values <- column_data[abs(modified_z_score) <= 3.5]
  
  # --- Continue as before, using only non-outlier values ---
  global_min <- min(non_outlier_values, na.rm = TRUE)
  if (global_min == Inf) global_min <- 1
  
  global_max <- max(non_outlier_values, na.rm = TRUE)
  if (global_max == -Inf) global_max <- 1
  
  ylim <- ceiling(global_max + abs(ceiling(0.1 * global_min)))
  ymin <- floor(global_min - abs(ceiling(0.1 * global_min)))
  abs_max <- max(abs(ymin), abs(ylim))
  
  abs_max <- if (abs_max >= 10) {
    round_10(abs_max)
  } else {
    round_5(abs_max)
  }
  
  # Determine step size and generate the yintercepts of the horizontal lines
  step_size <- case_when(
    abs_max * 2 <= 10 ~ 2,
    abs_max * 2 <= 25 ~ 5,
    abs_max * 2 <= 50 ~ 10,
    TRUE ~ 15
  )
  
  hlines <- sort(unique(c(seq(0, abs_max, by = step_size),
                          seq(0, abs_max, by = step_size) * -1)))
  
  list(hlines, abs_max)
}

# Compare two SummarizedExperiment conditions in a genomic window.
#
# Input:
#   - se_cond_1, se_cond_2: SummarizedExperiment objects with at least:
#       rowRanges() containing metadata columns:
#           - position
#           - half_life_dif
#           - log2_intensity_fold_change
#           - HL_dif_segment
#           - intensity_log2_fc_segment
#   - start_p: Integer, genomic start position (default: 0)
#
# Output:
#   - A list of ggplot2 objects for later PDF generation:
#       - diff_plot_hl_sense
#       - diff_plot_intensity_sense
#       - diff_plot_hl_antisense
#       - diff_plot_intensity_antisense


Cond_comp <- function(se_cond_1, se_cond_2, start_p = 0){
  end_p = start_p + 10000
  #Preparation of subsetting
  se_subset_sense <- GRanges(seqnames="chr",ranges = (start_p:end_p),  strand = "+")
  se_subset_antisense <- GRanges(seqnames="chr",ranges = (start_p:end_p),  strand = "-")
  
  ##################################################
  rR_sense <- se_cond_1 %>% subsetByOverlaps(se_subset_sense) %>% rowRanges()
  rR_antisense <- se_cond_1 %>% subsetByOverlaps(se_subset_antisense) %>% rowRanges()
  diff_df_sense <- subset(data.frame(half_life_dif = rR_sense$half_life_dif, log2_intensity_fold_change = rR_sense$log2_intensity_fold_change, position = rR_sense$position, HL_dif_segment = rR_sense$HL_dif_segment, intensity_log2_fc_segment = rR_sense$intensity_log2_fc_segment), position <= end_p)
  diff_df_antisense <- subset(data.frame(half_life_dif = rR_antisense$half_life_dif,log2_intensity_fold_change = rR_antisense$log2_intensity_fold_change, position = rR_antisense$position, HL_dif_segment = rR_antisense$HL_dif_segment, intensity_log2_fc_segment = rR_antisense$intensity_log2_fc_segment), position <= end_p)
  ##################################################
  # Half-life sense log2 fold change plot
  # Processing the segment and introducing possible outliers
  sense_hl_dif <- diff_df_sense %>%
    group_by(HL_dif_segment) %>%
    group_modify(~ process_segment(.x, "half_life_dif")) %>%
    mutate(
      min_pos = min(position, na.rm = TRUE),
      max_pos = max(position, na.rm = TRUE),
      mid_pos = min_pos + (max_pos - min_pos) / 2,
      mean_hl_dif = mean(half_life_dif, na.rm = TRUE)
    )
  # Get yintercepts of hlines
  hlines_info <- get_hlines(diff_df_sense, half_life_dif)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  # Cap outliers to within plotting limits
  sense_hl_dif <- sense_hl_dif %>%
    mutate(
      outliers = ifelse(!is.na(outliers) & abs(outliers) > abs_max,
                        sign(outliers) * abs_max,
                        outliers
      )
    )
  # Check if there are outliers in the data
  has_ol <- any(!is.na(sense_hl_dif$outliers))
  # Initialize plot object
  diff_plot_hl_sense <- ggplot(data = sense_hl_dif, aes(x = position, y = half_life_dif, color = as.factor(HL_dif_segment))) +
    geom_point(size=2.5) +
    {if(has_ol)geom_point(aes(x = position, y = outliers), shape = 13, size = 2.5)} +
    xlim(start_p,end_p) +
    geom_hline(yintercept=hlines, linetype = "dashed", color = "grey4") +
    geom_segment(data = sense_hl_dif, aes(x = min_pos, xend = max_pos, y=mean_hl_dif, yend = mean_hl_dif, color=as.factor(HL_dif_segment)),size = 0.7)+
    geom_text(
      data = distinct(sense_hl_dif, HL_dif_segment, .keep_all = TRUE),
      aes(x = mid_pos, y = mean_hl_dif, label = HL_dif_segment),
      color = "black",
      check_overlap = TRUE,
      family = "Roboto",
      size = 4,
      vjust = 3
    ) +
    scale_y_continuous(breaks = hlines, limits = c(-abs_max,abs_max)) +
    theme(legend.position="none",text = element_text(family = "Roboto"),panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank()) + labs(y = "half-life difference")
  # Intensity sense log2 fold change plot
  # Processing the segment and introducing possible outliers
  sense_int_fc <- diff_df_sense %>%
    group_by(intensity_log2_fc_segment) %>%
    group_modify(~ process_segment(.x, "log2_intensity_fold_change")) %>%
    mutate(
      min_pos = min(position, na.rm = TRUE),
      max_pos = max(position, na.rm = TRUE),
      mid_pos = min_pos + (max_pos - min_pos) / 2,
      mean_int_log2_fc = mean(log2_intensity_fold_change, na.rm = TRUE)
    )
  # Get yintercepts of hlines
  hlines_info <- get_hlines(sense_int_fc, log2_intensity_fold_change)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  sense_int_fc <- sense_int_fc %>%
    mutate(
      outliers = ifelse(!is.na(outliers) & abs(outliers) > abs_max,
                        sign(outliers) * abs_max,
                        outliers
      )
    )
  
  # Check if there are outliers in the data
  has_ol <- any(!is.na(sense_int_fc$outliers))
  # Initialize plot object
  diff_plot_int_sense <- ggplot(data = sense_int_fc, aes(x = position, y = log2_intensity_fold_change, color = as.factor(intensity_log2_fc_segment))) + geom_point(size=2.5) + xlim(start_p,end_p) +
    geom_point(size=2.5) +
    {if(has_ol)geom_point(aes(x = position, y = outliers), shape = 13, size = 2.5)} +
    xlim(start_p,end_p) +
    geom_segment(data = sense_int_fc, aes(x = min_pos, xend = max_pos, y=mean_int_log2_fc, yend = mean_int_log2_fc, color=as.factor(intensity_log2_fc_segment)),size = 0.7) +
    geom_hline(yintercept=hlines, linetype = "dashed", color = "grey4") +
    geom_text(
      data = distinct(sense_int_fc, intensity_log2_fc_segment, .keep_all = TRUE),
      aes(x = mid_pos, y = mean_int_log2_fc, label = intensity_log2_fc_segment),
      color = "black",
      check_overlap = TRUE,
      family = "Roboto",
      size = 4,
      vjust = 3
    ) +
    scale_y_continuous(breaks = hlines, limits = c(-abs_max,abs_max)) +
    theme(legend.position="none",text = element_text(family = "Roboto"),panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank()) + labs(y = "sense intensity log2 fc")
  
  # ANNOTATION-PLOT SENSE
  gff3_data <- gff3_preprocess(gff3_file_path)
  # Subsetting of gff3_data for sense strand and turning it into a "true" dataframe
  gff3_data_sense <- gff3_data %>%
    filter(strand == "+" & start <= end_p & end >= start_p) %>%
    mutate(
      start = pmax(start, start_p),  # Ensure start is not smaller than start_p
      end = pmin(end, end_p)         # Ensure end is not larger than end_p
    ) 
  # Formatting
  an_ylim = 10
  an_y_1 <- 0.25*an_ylim
  an_y_2 <- 0.75*an_ylim
  # Initialization of sense annotation plot object
  an_p_sense <- ggplot(gff3_data_sense, aes(x = start)) +
    ggchicklet:::geom_rrect(aes(xmin=start, xmax = end, ymin = an_y_1, ymax =an_y_2),color = "darkslategrey", fill = "azure3", alpha = 0.5, r = unit(0.5, 'npc')) + xlim(start_p,end_p) + ylim(0,an_ylim)+
    geom_text(data=gff3_data_sense, aes(x=start+(end-start)/2, y=0.5*an_ylim, label=gene_locus_modified),family = "Roboto", size=2.7,fontface = "bold") + labs(y = "+", face = "bold") +
    theme_void_yax()
  # theme_void_yax is also contained in the MB_plot_functions.R file!
  
  
  
  # half-life antisense plot
  # Processing the segment and introducing possible outliers
  antisense_hl_dif <- diff_df_antisense %>%
    group_by(HL_dif_segment) %>%
    group_modify(~ process_segment(.x, "half_life_dif")) %>%
    group_by(HL_dif_segment) %>%
    mutate(
      min_pos = min(position, na.rm = TRUE),
      max_pos = max(position, na.rm = TRUE),
      mid_pos = min_pos + (max_pos - min_pos) / 2,
      mean_hl_dif = mean(half_life_dif, na.rm = TRUE)
    )
  # Get yintercepts of hlines
  hlines_info <- get_hlines(diff_df_antisense, half_life_dif)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  antisense_hl_dif <- antisense_hl_dif %>%
    mutate(
      outliers = ifelse(!is.na(outliers) & abs(outliers) > abs_max,
                        sign(outliers) * abs_max,
                        outliers
      )
    )
  # Check if there are outliers in the data
  has_ol <- any(!is.na(antisense_hl_dif$outliers))
  # Initialize plot object
  diff_plot_hl_antisense <- ggplot(data = antisense_hl_dif, aes(x = position, y = half_life_dif, color = as.factor(HL_dif_segment))) +
    geom_point(size=2.5) +
    {if(has_ol)geom_point(aes(x = position, y = outliers), shape = 13, size = 2.5)} +
    geom_segment(data = antisense_hl_dif, aes(x = min_pos, xend = max_pos, y=mean_hl_dif, yend = mean_hl_dif, color=as.factor(HL_dif_segment)),size = 0.7) +
    xlim(start_p,end_p) +
    scale_y_continuous(breaks = hlines, limits = c(-abs_max,abs_max)) +
    geom_hline(yintercept=hlines, linetype = "dashed", color = "grey4") +
    geom_text(
      data = distinct(antisense_hl_dif, HL_dif_segment, .keep_all = TRUE),
      aes(x = mid_pos, y = mean_hl_dif, label = HL_dif_segment),
      color = "black",
      check_overlap = TRUE,
      family = "Roboto",
      size = 4,
      vjust = 3
    ) +
    theme(legend.position="none",text = element_text(family = "Roboto"),panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank()) + labs(y = "half-life difference")
  
  # Intensity antisense log2 fold change plot
  # Processing the segment and introducing possible outliers
  antisense_int_fc <- diff_df_antisense %>%
    group_by(intensity_log2_fc_segment) %>%
    group_modify(~ process_segment(.x, "log2_intensity_fold_change")) %>%
    mutate(
      min_pos = min(position, na.rm = TRUE),
      max_pos = max(position, na.rm = TRUE),
      mid_pos = min_pos + (max_pos - min_pos) / 2,
      mean_int_log2_fc = mean(log2_intensity_fold_change, na.rm = TRUE)
    )
  # Get yintercepts of hlines
  hlines_info <- get_hlines(antisense_int_fc, log2_intensity_fold_change)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  abs_calc = abs_max*0.35
  antisense_int_fc <- antisense_int_fc %>%
    mutate(
      outliers = ifelse(!is.na(outliers) & abs(outliers) > abs_max,
                        sign(outliers) * abs_max,
                        outliers
      )
    )
  # Check if there are outliers in the data
  has_ol <- any(!is.na(antisense_int_fc$outliers))
  # Initialize plot object
  diff_plot_int_antisense <- ggplot(data = antisense_int_fc, aes(x = position, y = log2_intensity_fold_change, color = as.factor(intensity_log2_fc_segment))) +
    geom_point(size=2.5) +
    {if(has_ol)geom_point(aes(x = position, y = outliers), shape = 13, size = 2.5)} +
    xlim(start_p,end_p) +
    geom_segment(data = antisense_int_fc, aes(x = min_pos, xend = max_pos, y=mean_int_log2_fc, yend = mean_int_log2_fc, color=as.factor(intensity_log2_fc_segment)),size = 0.7) +
    scale_y_continuous(breaks = hlines, limits = c(-abs_max,abs_max)) +
    geom_hline(yintercept=hlines, linetype = "dashed", color = "grey4") +
    geom_text(
      data = distinct(antisense_int_fc, intensity_log2_fc_segment, .keep_all = TRUE),
      aes(x = mid_pos, y = mean_int_log2_fc, label = intensity_log2_fc_segment),
      color = "black",
      check_overlap = TRUE,
      family = "Roboto",
      size = 4,
      vjust = 3
    ) +
    theme(legend.position="none",text = element_text(family = "Roboto"),panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank()) + labs(y = "antisense intensity log2 fc")
  
  # ANNOTATION PLOT ANTISENSE
  # Subsetting of gff3_data for sense strand and turning it into a "true" dataframe
  gff3_data_antisense <- gff3_data %>%
    filter(strand == "-" & start <= end_p & end >= start_p) %>%
    mutate(
      start = pmax(start, start_p),  # Ensure start is not smaller than start_p
      end = pmin(end, end_p)         # Ensure end is not larger than end_p
    ) 
  # Only Gene or Locus tag, not both
  # Initialization of antisense annotation plot object
  an_p_antisense <- ggplot(gff3_data_antisense, aes(x = start)) +
    ggchicklet:::geom_rrect(aes(xmin=start, xmax = end, ymin = an_y_1, ymax =an_y_2),color = "darkslategrey", fill = "azure3", alpha = 0.5, r = unit(0.5, 'npc')) + xlim(start_p,end_p) + ylim(0,an_ylim)+
    geom_text(data=gff3_data_antisense, aes(x=start+(end-start)/2, y=0.5*an_ylim, label=gene_locus_modified), size=2.7 ,family = "Roboto",fontface = "bold") + labs(y = "-",face = "bold") +
    theme_void_yax()
  
  
  
  print(ggarrange(diff_plot_int_sense, diff_plot_hl_sense,an_p_sense,an_p_antisense, diff_plot_hl_antisense,diff_plot_int_antisense, heights = c(2,2,1,1,2,2),ncol = 1, nrow = 6))
}

# ============================================
# final_vis_seg_comp
# Input: 
#   - se_cond_1: SummarizedExperiment for condition 1
#   - se_cond_2: SummarizedExperiment for condition 2
#   - name: filename for PDF output (default: "new_plot")
#   - gff3_file_path: path to GFF3 annotation file
#   - start_p: start position for plotting (default: 0)
#   - seg_end: end position (default: max(se_cond_1 positions))
# Output: 
#   - Saves multi-page PDF with condition comparison plots
# ============================================

final_vis_seg_comp <- function(se_cond_1, se_cond_2, name = "new_plot",gff3_file_path, start_p = 0 , seg_end = max(rowRanges(se_cond_1)$position, na.rm = TRUE)) {
  on.exit({
    while (dev.cur() > 1) dev.off()
    final_message <- "You are finally done after an eternity of waiting."
    print(final_message)
  }, add = TRUE)
  start_time <- Sys.time()
  # Initialize pdf
  pdf(paste0("~/Desktop/R_plots/", name, ".pdf"), width = 16, height = 10)
  
  # Function to generate plot for each page
  generate_plot_page_comp <- function(start_p) {
    suppressMessages({
      Cond_comp(se_cond_1,se_cond_2, start_p)
      
      ind_seg_end <- start_p + 10000
      progress_message <- paste0(
        "Plotting position ", start_p, " to position ", ind_seg_end,
        " (Page ", (ind_seg_end / 10000), "). Roughly at ",
        ceiling((start_p / seg_end) * 100), "% completion."
      )
      print(progress_message)
    })
  }
  
  # Loop through segments
  positions <- seq(start_p, seg_end, by = 10000)
  lapply(positions, generate_plot_page_comp)
  
  # Close PDF device and print the time it took for the function to run
  dev.off()
  end_time <- Sys.time()
  elapsed_t <- end_time - start_time
  print(elapsed_t)
}

