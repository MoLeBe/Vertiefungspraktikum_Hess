library("SummarizedExperiment")
library("RcppDynProg")
library("GenomicRanges")
library("ggplot2")
library("ggridges")
library("dplyr")
library("foreach")
library("doParallel")
library("SummarizedExperiment")
library("parallel")
library("stringr")
library("rtracklayer")
library("ggpubr")
library("scales")
library("ggchicklet")
library("ggrepel")
library("ggsci")
library("rlang")
library("readr")
library("data.table")
# Add Roboto font (download if needed)
font_add_google("Roboto", "Roboto")

# Enable showtext
showtext_auto()

# ============================================
# seg_pen
# Input: 
#   - inp: data frame with RNA-seq data
#   - test_type: character, type of statistical test ("wx" for Wilcoxon, "t" for t-test)
#   - cores: number of cores for parallel processing
#   - col_name: column name for analysis (e.g., "half_life" or "intensity")
#   - penalty_vector: vector of penalty values for segmentation
#   - beta: constant used in F-value calculation (default = 5)
# Output: list with:
#   - segmentation_result_dataframe: data frame of segment results
#   - penalty_quality_dataframe: data frame with precision, recall, and F-value for each penalty
# ============================================

# PENALTY DETERMINATION
# Main function to perform the segmentation and analysis
seg_pen <- function(inp, test_type = "wx", cores = detectCores() - 2, col_name = "half_life",  penalty_vector , beta = 5){
  start_time <- Sys.time() 
  # Register parallel backend to use multiple cores
  registerDoParallel(cores)
  registerDoMC(cores)
  
  # Sort the dataframe by strand and position
  inp <- rifi:::inp_order(inp)
  
  # Create a temporary dataframe with delay segments
  tmp_df <- inp_df(inp, "ID", col_name, "position_segment")
  
  # Reverse the order for the minus strand
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  # Clean up the position_segment column
  tmp_df[, "position_segment"] <- gsub("_O|_NA", "", tmp_df$position_segment)
  
  # Remove rows with NA values
  tmp_df <- na.omit(tmp_df)
  
  # Get unique segments
  unique_seg <- unique(tmp_df$position_segment)
  
  # Dynamic Programming: Prepare segments
  segs <- mclapply(unique(tmp_df[["position_segment"]]), function(seg) {
    tmp_df[tmp_df[["position_segment"]] == seg, ]
  }, mc.cores = cores)
  # Function to process segments with dynamic programming and Wilcoxon tests
  penalty_seg_processing <- function(segment, pen) {
    # Initialize a segment counter
    number <- 0
    
    # Convert the segment into a dataframe
    tmp_df <- data.frame(segment)
    
    # Generate random half-life values (for debugging purposes)
    #tmp_df$half_life <- runif(length(tmp_df$half_life), 1, 50)
    
    # Use dynamic programming to find segment limits using an Rcpp function
    col_data <- tmp_df[[col_name]]
    segment_tmp_df <- solve_for_partition(1:nrow(tmp_df), col_data, penalty = pen)
    
    # Extract segment limits
    seg_range_list <- segment_tmp_df[["x"]]
    seg_seq <- seq(1, length(seg_range_list), 2)
    
    # Process each segment
    segment_l <- mclapply(seg_seq, function(k) {
      # Get the range for the current segment
      if (length(seg_range_list) > 1) {
        segment_range <- seg_range_list[k]:seg_range_list[k + 1]
      } else {
        segment_range <- seg_range_list[k]
      }
      
      # Increment the segment counter
      number <- number + 1
      
      # Subset the temporary dataframe to get the data of the current segment
      segment_df <- tmp_df[segment_range, ]
      
      # Add columns to the dataframe with the segment identifier and penalty
      segment_df[["position_segment"]] <- paste0("Dc_", number + tmp_df[["position"]][1])
      segment_df[["penalty"]] <- pen
      segment_df[["test"]] <- NA
      segment_df[["p_value"]] <- NA
      
      return(segment_df)
    }, mc.cores = cores)
    # Perform pairwise Wilcoxon or t-tests if there is more than one segment
    if (length(segment_l) > 1) {
      for (i in 2:length(segment_l)) {
        test_data <- list(segment_l[[i]][[col_name]], segment_l[[i - 1]][[col_name]])
        test_results <- switch(test_type,
                               "wx" = wilcox.test(test_data[[1]], test_data[[2]], alternative = "two.sided"),
                               "t"  = t.test(test_data[[1]], test_data[[2]], alternative = "two.sided", var.equal = FALSE)
        )
        segment_l[[i]][["test"]][1] <- test_results[["p.value"]] < 0.05
        segment_l[[i]][["p_value"]][1] <- test_results[["p.value"]]
      }
    }
    # Combine all segment results into a single dataframe
    combined_segments <- bind_rows(segment_l)
    return(combined_segments)
  }
  
  
  # Apply the penalty_seg_processing function for each penalty
  uni_segs <- mclapply(penalty_vector, function(pen) {
    segs_p <- mclapply(segs, penalty_seg_processing, pen, mc.cores = cores)
    segs_p_combined <- do.call(rbind, segs_p)
    test_table <- table(segs_p_combined[["test"]], useNA = "ifany")
    segs_p_combined <- segs_p_combined %>% mutate(correct_segments = as.integer(test_table["TRUE"]),
                          incorrect_segments = as.integer(test_table["FALSE"]))
    return(segs_p_combined)
  }, mc.cores = cores)                 
  
  
  uni_results_df <- bind_rows(uni_segs)
  pen_quality_df <- uni_results_df %>%
    group_by(penalty) %>%
    summarize(correct_segments = first(correct_segments),
              incorrect_segments = first(incorrect_segments),
              .groups = "drop") %>%
    mutate(Precision = correct_segments / (correct_segments + incorrect_segments),
           Recall = correct_segments / max(correct_segments),
           F_value = 2 * (((1 + beta^2) * Precision * Recall) / ((beta^2 * Precision) + Recall)))
  
  uni_results_df <- as.data.frame(bind_rows(uni_segs))
  end_time <- Sys.time()
  elapsed_t <- end_time - start_time
  print(elapsed_t) 
  # Return the final data frame as the function's output
  return(list(
    segmentation_result_dataframe = uni_results_df,
    penalty_quality_dataframe = pen_quality_df
  ))
}

# ============================================
# fine_tune_penalty
# Input: 
#   - inp: data frame with RNA-seq data
#   - col_name: column name for analysis (e.g., "half_life" or "intensity")
#   - initial_penalty: initial penalty value for fine-tuning
# Output: list with:
#   - penalty_quality_dataframe: data frame with precision, recall, and F-value for each penalty
#   - best_penalty: the penalty value with the highest F-value
# ============================================

fine_tune_penalty <- function(inp, col_name, initial_penalty) {
  fine_penalty_vector <- seq(initial_penalty - 5, initial_penalty + 5)
  pen_output_fine <- seg_pen(inp, test_type = "wx", col_name = col_name, penalty_vector = fine_penalty_vector)
  pen_quality_df_fine <- pen_output_fine[["penalty_quality_dataframe"]]
  best_penalty <- pen_quality_df_fine %>% filter(F_value == max(F_value)) %>% pull(1) %>% { if (length(.) > 1) .[1] else .}
  return(list(penalty_quality_dataframe = pen_quality_df_fine, best_penalty = best_penalty))
}
# ============================================
# penalty_det_wrapper
# Input: 
#   - inp: data frame with RNA-seq data
#   - test_type: character, type of statistical test ("wx" for Wilcoxon, "t" for t-test)
#   - penalty_vector: vector of penalty values for segmentation
# Output: list with:
#   - full_penalty_list: list with best penalties for half-life and intensity
#   - hl_pen_quality_df_broad: data frame with precision, recall, and F-value for half-life penalties
#   - hl_pen_quality_df_fine: data frame with fine-tuned precision, recall, and F-value for half-life penalties
#   - hl_best_penalty: best penalty for half-life segmentation
#   - int_pen_quality_df_broad: data frame with precision, recall, and F-value for intensity penalties
#   - int_pen_quality_df_fine: data frame with fine-tuned precision, recall, and F-value for intensity penalties
#   - int_best_penalty: best penalty for intensity segmentation
# ============================================

penalty_det_wrapper <- function(inp, test_type = "wx", penalty_vector = seq(20,25,5)){
  # Determine initial and fine-tuned half-life penalties
  hl_pen_output_broad <- seg_pen(inp, test_type = test_type, penalty_vector = penalty_vector)
  hl_pen_quality_df_broad <- hl_pen_output_broad[["penalty_quality_dataframe"]]
  hl_initial_penalty <- hl_pen_quality_df_broad %>% filter(F_value == max(F_value)) %>% pull(1) %>% { if (length(.) > 1) .[1] else .}
  hl_pen_output_fine <- fine_tune_penalty(inp, "half_life", hl_initial_penalty)
  hl_pen_quality_df_fine <- hl_pen_output_fine[["penalty_quality_dataframe"]]
  hl_best_penalty <- hl_pen_output_fine[["best_penalty"]]
  # Determine initial and fine-tuned intensity penalties
  int_pen_output_broad <- seg_pen(inp, test_type = test_type, col_name = "intensity", penalty_vector = penalty_vector)
  int_pen_quality_df_broad <- int_pen_output_broad[["penalty_quality_dataframe"]]
  int_initial_penalty <- int_pen_quality_df_broad %>% filter(F_value == max(F_value)) %>% pull(1) %>% { if (length(.) > 1) .[1] else .}
  int_pen_output_fine <- fine_tune_penalty(inp, "intensity", int_initial_penalty)
  int_pen_quality_df_fine <- int_pen_output_fine[["penalty_quality_dataframe"]]
  int_best_penalty <- int_pen_output_fine[["best_penalty"]]
  # Combine the best penalties
  full_penalty_list <- list(best_half_life_penalty = hl_best_penalty, best_intensity_penalty = int_best_penalty)
  return(list(full_penalty_list = full_penalty_list, 
    hl_pen_quality_df_broad = hl_pen_quality_df_broad,
    hl_pen_quality_df_fine = hl_pen_quality_df_fine,
    hl_best_penalty = hl_best_penalty,
    int_pen_quality_df_broad = int_pen_quality_df_broad,
    int_pen_quality_df_fine = int_pen_quality_df_fine,
    int_best_penalty = int_best_penalty
  ))
}




############################################################################################################
# SEGMENTATION
# ============================================
# Rcpp_segment
# Input: 
#   - inp: SummarizedExperiment object
#   - col_name: Character string, one of "half_life" or "intensity"
#   - pen: Integer, penalty value used for segmentation
# Output: 
#   - inp_filtered: SummarizedExperiment object with updated segment identifiers for half life or intensity
# ============================================

Rcpp_segment <- function(inp, col_name = "half_life", pen = 20) {
  start_time <- Sys.time() 
  # I.Preparations: the dataframe is configured and some other variables are
  seg_name <- switch(col_name,
                     "half_life" = "HL_fragment",
                     "intensity" = "intensity_fragment"
  )
  identifier <- switch(col_name,
                       "half_life" = "Dc_",
                       "intensity" = "I_"
  )
  base_seg <- switch(col_name,
                     "half_life" = "position_segment",
                     "intensity" = "HL_fragment"
  )
  mcols(rowRanges(inp))[[seg_name]] <- NA
  # the dataframe is sorted by strand and position.
  inp <- rifi:::inp_order(inp) 
  
  tmp_df <- inp_df(inp, "ID", col_name, base_seg)
  # Log-transform the intensity data if applicable
  if (col_name == "intensity") {
    tmp_df[[col_name]] <- log2(tmp_df[[col_name]] + 1)  # log2 transformation; add 1 to avoid log(0)
  }
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  unique_seg <- unlist(unique(tmp_df[[base_seg]]))
  
  count <- 1
  # II. Dynamic Programming: the scoring function is interpreted
  
  segs <- mclapply(unique(tmp_df[[base_seg]]), function(seg) {
    tmp_df[tmp_df[[base_seg]] == seg, ]
  }, mc.cores = cores)
  
  
  Rcpp_seg_df <- data.frame(ID = numeric(),
                            FLT = numeric(), 
                            strand = character(), 
                            position = numeric())
  Rcpp_seg_df[[base_seg]] <- character()
  Rcpp_seg_df[[col_name]] <- numeric()
  
  number = 0
  
  Rcpp_process <- function(segment, pen){
    # Temporary dataframe containing the data of the current segment
    tmp_df <- data.frame(segment)
    # Temporary dataframe containing the Rcpp segmented data
    Rcpp_tmp_df <- solve_for_partition(1:nrow(tmp_df), tmp_df[[col_name]], penalty = pen,min_seg = 3)
    
    
    # Extract the segment limits
    seg_list <- Rcpp_tmp_df[["x"]]
    # Generate sequence indices
    seg_seq <- seq(from = 1, to = length(seg_list), by = 2)
    result <- lapply(seg_seq, function(k) {
      if (length(seg_list) != 1) {
        # Increment the segment counter
        number <<- number + 1
        # Get the range for the current segment
        segment_range <- seg_list[k]:seg_list[k + 1]
      } else {
        # Increment the segment counter
        number <<- number + 1
        # Get the range for a single-element segment
        segment_range <- seg_list[k]
      }
      
      # Subset the temporary dataframe to get the data of the current segment
      Rcpp_tmp_df_subset <- tmp_df[segment_range, ]
      # Add a column with the segment identifier
      Rcpp_tmp_df_subset[["segment_ID"]] <- switch(col_name,
                                              "half_life" = paste0("Dc_", number),
                                              "intensity" = paste0("I_", number)
      )
      
      # Return the subset data frame for the current segment
      Rcpp_tmp_df_subset
    })
    do.call(rbind, result)
  }
  final_segs <- lapply(segs, Rcpp_process, pen)
  final_segs_df <- do.call(rbind, final_segs)
  matching_rows <- rownames(inp) %in% rownames(final_segs_df)
  # Subset the inp object based on matching rows
  inp_filtered <- inp[matching_rows, ]
  mcols(rowRanges(inp_filtered))[[seg_name]] <- final_segs_df$segment_ID
  end_time <- Sys.time() 
  elapsed_t <- end_time - start_time
  print(elapsed_t)
  inp_filtered
}


############################################################################################################
# VISUALIZATION
# Colors
locuszoom_colors <- pal_locuszoom("default")(7)
extended_locuszoom_colors <- rep(locuszoom_colors, length.out = 100)
# Define own annotation plot theme
theme_void_yax <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(face = "bold",family = "Roboto"),
        legend.position="none",
        panel.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
}
# ============================================
# gff3_preprocess
# Input: 
#   - path: Character string, file path to a GFF3 file
#   - next_line: Logical, whether to replace underscores in gene names with line breaks
# Output: 
#   - tmp: Dataframe with processed GFF3 data, including gene locus modifications
# ============================================


# Processing function for gff3 file to get it in appropriate format
gff3_preprocess <- function(path, next_line = TRUE) {
  
  # Read the GFF file and select necessary columns
  inp <- readGFF(path)
  tmp <- as.data.frame(inp[, c("type", "start", "end", "strand", "gene", "locus_tag")])
  
  # Ensure columns are character vectors
  tmp[["gene"]] <- as.character(tmp[["gene"]])
  tmp[["locus_tag"]] <- as.character(tmp[["locus_tag"]])
  
  # Remove quotes from gene and locus_tag columns
  tmp[["gene"]] <- gsub('"', "", tmp[["gene"]], fixed = TRUE)
  tmp[["locus_tag"]] <- gsub('"', "", tmp[["locus_tag"]], fixed = TRUE)
  
  # Filter rows based on type
  tmp <- tmp %>%
    filter(str_detect(type, "^CDS$|UTR|asRNA|antisense_RNA|ncRNA|^tRNA$"))
  
  colnames(tmp)[1] <- "region"
  
  # Replace NA or "NA" in gene with locus_tag
  tmp <- tmp %>%
    mutate(gene_locus = ifelse(is.na(gene) | gene == "NA", locus_tag, gene)) %>%
    mutate(
      gene_locus_modified = case_when(
        next_line == TRUE ~ case_when(
          str_detect(gene_locus, "_") ~ str_replace(gene_locus, "_", "\n"),
          str_length(gene_locus) > 5 & str_detect(gene_locus, "\\d") ~ str_replace(gene_locus, "(\\d)", "\n\\1"),
          TRUE ~ gene_locus
        ),
        TRUE ~ gene_locus  # Default case if next_line is not TRUE
      ),
      gene_locus_modified = gene_locus_modified %>%
        str_replace_all("^c\\(|\\)$", "") %>%
        str_replace_all("\\s*,\\s*", ",") %>%
        str_replace_all("\\\\", "")
    )
  return(tmp)
}
# ============================================
# get_hlines
# Input: 
#   - df: Dataframe with plot data
#   - y_col: Column name for which the horizontal lines should be calculated
# Output: 
#   - hlines: Vector of y-intercepts for horizontal lines
#   - abs_max: Maximum y-value used for axis limits
# ============================================
# Function for definition of yintercepts of hlines in the plot
get_hlines <- function(df,y_col){
  y_col <- as.character(substitute(y_col))
  column_data <- df %>% pull(!!sym(y_col))
  global_min <- min(column_data, na.rm = TRUE)
  {if(global_min==Inf)global_min = 1}
  global_max <- max(column_data,na.rm = TRUE)
  {if(global_max==-Inf)global_max = 1}
  ylim <- ceiling(global_max + abs(ceiling(0.1 * global_min)))
  ymin <- floor(global_min - abs(ceiling(0.1 * global_min)))
  abs_max <- max(abs(ymin), abs(ylim))
  x <- if(abs_max < 20){
    round_5(abs_max)
  } else{
    abs_max = abs_max
  }
  
  # Determine step size and generate the yintercepts of the horizontal lines
  step_size <- case_when(
    abs_max*2 <= 10 ~ 2,
    abs_max*2 <= 42 ~ 5,
    abs_max*2 <= 60 ~ 10,
    TRUE ~ 15
  )
  hlines <- sort(unique(c(seq(0,abs_max, by = step_size),seq(0,abs_max, by = step_size)*-1)))
  list(hlines, abs_max)
}
# ===========================================
# process_segment_Rcpp
# Input: 
#   - segment_data: Dataframe with segment data
#   - column_name: Character string, column name to process (e.g., "half_life")
# Output: 
#   - segment_data: Dataframe with outliers identified and replaced with NA
# ===========================================
process_segment_Rcpp <- function(segment_data, column_name) {
  
  # Pull the specified column from the data
  column_data <- segment_data[[column_name]]
  
  # Calculate the median and MAD (median absolute deviation)
  mad_value <- mad(column_data, constant = 1.4826, na.rm = TRUE)
  
  # Ensure mad_value is not NA or 0
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- 1
  }
  
  median_value <- median(column_data, na.rm = TRUE)
  
  # Apply the Modified Z-score method to identify outliers
  segment_data <- segment_data %>%
    mutate(
      modified_z_score = 0.6745 * (column_data - median_value) / mad_value,
      outliers = ifelse(abs(modified_z_score) > 3.5, !!sym(column_name), NA),
      !!sym(column_name) := ifelse(abs(modified_z_score) > 3.5, NA, !!sym(column_name))
    ) %>%
    select(-modified_z_score) # Remove the modified_z_score column if no longer needed
  
  return(segment_data)
}
# ===========================================
# segmentation_vis_parts
# Input: 
#   - inp_se: SummarizedExperiment object
#   - start_p: Integer, starting position for the plot
# Output: 
#   - hl_p_sense: ggplot object for half-life plot (sense strand)
#   - int_p_sense: ggplot object for intensity plot (sense strand)
#   - an_p_sense: ggplot object for annotation plot (sense strand)
#   - hl_p_antisense: ggplot object for half-life plot (antisense strand)
#   - int_p_antisense: ggplot object for intensity plot (antisense strand)
#   - an_p_antisense: ggplot object for annotation plot (antisense strand)
# These are used in the final_vis_seg function to generate the pdf
# ===========================================
segmentation_vis_parts <- function(inp_se, start_p = 0){
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
    group_by(HL_segment) %>%
    group_modify(~ process_segment_Rcpp(.x, "half_life")) %>%
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
    group_by(HL_segment) %>%
    group_modify(~ process_segment_Rcpp(.x, "half_life")) %>%
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

# ===========================================
# final_vis_seg
# Input:
#   - inp: A SummarizedExperiment object containing the data to be visualized.
#   - name: A string specifying the name for the output PDF file. Default is "new_plot".
#   - seg_end: An integer specifying the end position for segmentation. Default is the maximum position from rowRanges(inp).

# Output:
#   - A PDF file containing multiple pages, each with a plot segment of 10000 nt.
#   - The PDF is saved on the user's desktop with the name specified by name (default is "new_plot.pdf").
# ===========================================

final_vis_seg <- function(inp, name = "new_plot", seg_end = max(rowRanges(inp)$position, na.rm = TRUE)) {
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
      segmentation_vis_parts(inp, start_p)
      
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
