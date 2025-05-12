
# ============================================
# final_vis_seg_rifid
# Input: Summarized Experiment object containing rifi segmentation data
# Output: 
#   - Saves multi-page PDF with rifi segmentation plots
# ============================================
library(dplyr)
library(rlang)
library(stringr)
# Add Roboto font (download if needed)
font_add_google("Roboto", "Roboto")

# Enable showtext
showtext_auto()

process_segment_rifi <- function(segment_data, column_name) {
  
  # Determine the identifier column based on the column_name
  identifier_col <- "HL_segment"
  
  segment_data <- segment_data %>%
    mutate(
      # Check if identifier ends in "_O"
      is_outlier = str_ends(.data[[identifier_col]], "_O"),
      # Remove "_O" from identifier if it’s there
      !!identifier_col := ifelse(is_outlier,
                                 str_remove(.data[[identifier_col]], "_O$"),
                                 .data[[identifier_col]]),
      # Store the outlier value in a new column, otherwise NA
      outliers = ifelse(is_outlier, .data[[column_name]], NA),
      # Set the original column to NA for outliers
      !!sym(column_name) := ifelse(is_outlier, NA, .data[[column_name]])
    ) %>%
    select(-is_outlier)
  
  return(segment_data)
}




segmentation_vis_parts_rifi <- function(inp_se, start_p = 0){
  end_p = start_p + 10000
  # Subsetting of summarized Experiment object to desired length and strand
  se_subset_sense <- GRanges(
    seqnames = "chr",
    ranges = IRanges(start = start_p, end = end_p),
    strand = "+"
  )
  se_subset_antisense <- GRanges(
    seqnames = "chr",
    ranges = IRanges(start = start_p, end = end_p),
    strand = "-"
  )
  rR_sense <- inp_se %>% subsetByOverlaps(se_subset_sense) %>% rowRanges()
  rR_antisense <- inp_se %>% subsetByOverlaps(se_subset_antisense) %>% rowRanges()
  
  # HALF-LIFE SENSE
  # Initialization of sense dataframe and subsetting to appropriate length
  ggplot_df_sense <- subset(data.frame(half_life = rR_sense$half_life,
                                       intensity = rR_sense$intensity,
                                       position = rR_sense$position,
                                       HL_segment = rR_sense$HL_fragment,
                                       intensity_segment = rR_sense$intensity_fragment),
                            position <= end_p)
  
  
  # Calculate the global maximum of the half_life and the ylim_hl_sense
  global_max_hl <- ceiling(max(ggplot_df_sense$half_life, na.rm = TRUE))
  ylim_hl_sense <- global_max_hl + ceiling(0.1 * global_max_hl)
  
  # Apply the function to each segment
  ggplot_df_sense <- ggplot_df_sense %>%
    process_segment_rifi("half_life")
  ggplot_df_sense <- ggplot_df_sense %>%
    group_by(HL_segment) %>%
    mutate(mean_hl = mean(half_life, na.rm = TRUE),
           min_pos_hl = min(position, na.rm = TRUE),
           max_pos_hl = max(position, na.rm = TRUE),
           mid_pos_hl = min_pos_hl + (max_pos_hl - min_pos_hl) / 2) %>%
    ungroup() %>%
    group_by(intensity_segment)  %>%
    mutate(mean_int = mean(intensity, na.rm = TRUE),
           min_pos_int = min(position, na.rm = TRUE),
           max_pos_int = max(position, na.rm = TRUE),
           mid_pos_int = min_pos_int + (max_pos_int - min_pos_int) / 2) %>%
    ungroup()
  ggplot_df_sense <- ggplot_df_sense %>%
    mutate(
      half_life = pmin(half_life, 20),
      outliers = if ("outliers" %in% names(.)) pmin(outliers, 20) else NA
    )
  # Determine step size and generate the yintercepts of the horizontal lines
  hlines_info <- get_hlines(ggplot_df_sense, half_life)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  # Determine if there are outliers in the data
  has_ol <- any(!is.na(ggplot_df_sense$outliers))
  # Initialization of plot object
  hl_p_sense <- ggplot(data=ggplot_df_sense, aes(x = position,
                                                 y = half_life,
                                                 color = as.factor(HL_segment))) +
    geom_segment(data = ggplot_df_sense, aes(x = min_pos_hl, xend = max_pos_hl, y=mean_hl, yend = mean_hl, color=as.factor(HL_segment)),size = 0.7) +
    geom_point(size=2.5) + 
    {if(has_ol)geom_point(aes(x = position, y = outliers),
                          shape = 13, size = 2.5)} + 
    xlim(start_p,end_p) + 
    scale_y_continuous(breaks = hlines,
                       limits = c(0,abs_max)) + 
    geom_hline(yintercept=hlines, 
               linetype = "dashed", color = "grey4") +
    scale_color_manual(values = extended_locuszoom_colors) +
    theme(legend.position="none",
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) +
    labs(y = "half-life [min]", family = "Roboto")
  # INTENSITY-LIFE SENSE
  min_y_int_sense <- floor(min(ggplot_df_sense$intensity)) # Determine lowest intensity, that has to be plotted
  hlines_int_sense <- c(2^7,2^9,2^11,2^13,2^15)
  # Initialization of plot object
  int_p_sense <- ggplot(data=ggplot_df_sense, aes(x = position,
                                                  y = intensity, 
                                                  color = as.factor(intensity_segment))) +
    geom_segment(data = ggplot_df_sense, aes(x = min_pos_int, xend = max_pos_int, y=mean_int, yend = mean_int, color=as.factor(intensity_segment)),size = 0.7) +
    geom_point(size=2.5) + 
    xlim(start_p,end_p)  + 
    ylim(2^7,2^15) +
    geom_hline(yintercept=hlines_int_sense,
               linetype = "dashed", 
               color = "grey4") +
    scale_color_manual(values = extended_locuszoom_colors) +
    theme(legend.position="none",
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + 
    scale_y_continuous(trans='log2', 
                       breaks = hlines_int_sense,
                       labels = c(expression("2"^7),
                                  expression("2"^9),
                                  expression("2"^11),
                                  expression("2"^13), 
                                  expression("2"^15)))  + 
    labs(y = "intensity [A.U]")
  
  
  # ANNOTATION-PLOT SENSE
  # Read gff3_data for plotting and annotation of genes
  gff3_data <- gff3_preprocess(gff3_file_path)
  #Subsetting of gff3_data for sense strand and turning it into a "true" dataframe
  gff3_data_sense <- as.data.frame(subset(gff3_data, strand == "+" & start > start_p & end < end_p))
  # Annotation plot "architecture"
  an_ylim = 10
  an_y_1 <- 0.25*an_ylim 
  an_y_2 <- 0.75*an_ylim
  
  # Initialization of annotation plot object
  an_p_sense <- ggplot(gff3_data_sense, aes(x = start)) +
    ggchicklet:::geom_rrect(aes(xmin=start, xmax = end, ymin = an_y_1, ymax =an_y_2),color = "darkslategrey", fill = "azure3", alpha = 0.5, r = unit(0.5, 'npc')) + xlim(start_p,end_p) + ylim(0,an_ylim)+
    geom_text(data=gff3_data_sense, aes(x=start+(end-start)/2, y=0.5*an_ylim, label=gene_locus_modified), family = "Roboto", size=2.7,fontface = "bold") + labs(y = "+", face = "bold") +
    theme_void_yax()
  
  
  # HALF-LIFE ANTISENSE
  # Initialization of antisense dataframe and subsetting to appropriate length
  ggplot_df_antisense <- subset(data.frame(half_life = rR_antisense$half_life, intensity = rR_antisense$intensity, position = rR_antisense$position, HL_segment = rR_antisense$HL_fragment, intensity_segment = rR_antisense$intensity_fragment), position <= end_p)
  
  # Apply the function to each segment
  ggplot_df_antisense <- ggplot_df_antisense %>%
    process_segment_rifi("half_life")
  ggplot_df_antisense <- ggplot_df_antisense %>%
    group_by(HL_segment) %>%
    mutate(mean_hl = mean(half_life, na.rm = TRUE),
           min_pos_hl = min(position, na.rm = TRUE),
           max_pos_hl = max(position, na.rm = TRUE),
           mid_pos_hl = min_pos_hl + (max_pos_hl - min_pos_hl) / 2) %>%
    ungroup() %>%
    group_by(intensity_segment) %>%
    mutate(mean_int = mean(intensity, na.rm = TRUE),
           min_pos_int = min(position, na.rm = TRUE),
           max_pos_int = max(position, na.rm = TRUE),
           mid_pos_int = min_pos_int + (max_pos_int - min_pos_int) / 2) %>%
    ungroup()
  ggplot_df_antisense <- ggplot_df_antisense %>%
    mutate(
      half_life = pmin(half_life, 20),
      outliers = if ("outliers" %in% names(.)) pmin(outliers, 20) else NA
    )
  # Determine step size and generate the yintercepts of the horizontal lines
  hlines_info <- get_hlines(ggplot_df_antisense, half_life)
  hlines <- hlines_info[[1]]
  abs_max <- hlines_info[[2]]
  # Determine if there are outliers in the data
  has_ol <- any(!is.na(ggplot_df_antisense$outliers))
  # Initialization of plot object
  hl_p_antisense <- ggplot(data=ggplot_df_antisense, aes(x = position, y = half_life, color = as.factor(HL_segment))) + geom_point(size=2.5) + {if(has_ol)geom_point(aes(x = position, y = outliers), shape = 13, size = 2.5)} + xlim(start_p,end_p) + scale_y_reverse(breaks = hlines, limits = c(abs_max,0)) + geom_hline(yintercept=hlines, linetype = "dashed", color = "grey4") +
    geom_segment(data = ggplot_df_antisense, aes(x = min_pos_hl, xend = max_pos_hl, y=mean_hl, yend = mean_hl, color=as.factor(HL_segment)),size = 0.7) +
    scale_color_manual(values = extended_locuszoom_colors) +
    theme(legend.position="none",panel.background=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + labs(y = "half-life [min]") 
  # INTENSITY PLOT ANTISENSE
  min_y_int_antisense <- floor(min(ggplot_df_antisense$intensity)) 
  hlines_int_antisense <- c(2^7,2^9,2^11,2^13,2^15)
  # Define transformation function of y-axis (reverse because it is plotted upside down)
  trans_log2_rev <- trans_new("log2_reverse", transform = function(x) -log2(x),inverse = function(x) 2^(-x))
  # Initialization of plot object
  int_p_antisense <- ggplot(data=ggplot_df_antisense, aes(x = position, y = intensity, color = as.factor(intensity_segment))) + geom_point(size=2.5) + xlim(start_p,end_p) + ylim(2^7,2^15) + geom_hline(yintercept=hlines_int_antisense, linetype = "dashed", color = "grey4") +
    scale_color_manual(values = extended_locuszoom_colors) +
    geom_segment(data = ggplot_df_antisense, aes(x = min_pos_int, xend = max_pos_int, y=mean_int, yend = mean_int, color=as.factor(intensity_segment)),size = 0.7) +
    theme(legend.position="none",panel.background=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + labs(y = "intensity [A.U]") + scale_y_continuous(trans= trans_log2_rev, breaks = hlines_int_antisense, labels = c(expression("2"^7), expression("2"^9), expression("2"^11), expression("2"^13), expression("2"^15)))
  
  # ANNOTATION PLOT ANTISENSE
  # Subsetting of gff3_data for sense strand and turning it into a "true" dataframe
  gff3_data_antisense <- as.data.frame(subset(gff3_data, strand == "-" & start > start_p & end < end_p))
  # Initialization of annotation plot object
  an_p_antisense <- ggplot(gff3_data_antisense, aes(x = start)) +
    ggchicklet:::geom_rrect(aes(xmin=start, xmax = end, ymin = an_y_1, ymax =an_y_2),color = "darkslategrey", fill = "azure3", alpha = 0.5, r = unit(0.5, 'npc')) + xlim(start_p,end_p) + ylim(0,an_ylim)+
    geom_text(data=gff3_data_antisense, aes(x=start+(end-start)/2, y=0.5*an_ylim, label=gene_locus_modified),family = "Roboto", size=2.7,fontface = "bold") + labs(y = "-",face = "bold") +
    theme_void_yax()
  
  
  # Compilation of all plot objects initialized in the previous steps and printing to a page in the pdf-file!
  print(ggarrange(int_p_sense,hl_p_sense,an_p_sense,an_p_antisense,hl_p_antisense,int_p_antisense, heights = c(2,2,1,1,2,2),ncol = 1, nrow = 6))
  
}

final_vis_seg_rifi <- function(inp, name = "new_plot", seg_end = max(rowRanges(inp)$position, na.rm = TRUE)) {
  on.exit({
    while (dev.cur() > 1) dev.off()
    final_message <- "You are finally done after an eternity of waiting."
    print(final_message)
  }, add = TRUE)
  start_time <- Sys.time() 
  # Initialize pdf
  pdf(paste0("~/Desktop/R_plots/", name, ".pdf"), width = 16, height = 10)
  
  # Function to generate plot for each page
  generate_plot_page <- function(start_p) {
    suppressMessages({
      segmentation_vis_parts_rifi(inp, start_p)
      
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
  positions <- seq(0, seg_end, by = 10000)
  lapply(positions, generate_plot_page)
  
  # Close PDF device and print the time it took for the function to run
  dev.off()
  end_time <- Sys.time() 
  elapsed_t <- end_time - start_time
  print(elapsed_t)
}