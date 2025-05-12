# ============================================
# segment_HL_Rcpp_comp
# Input: 
#   - se_cond_1: SummarizedExperiment for condition 1
#   - se_cond_2: SummarizedExperiment for condition 2
#   - cores: number of cores for parallel processing (default: 1)
#   - pen: penalty value for dynamic programming (default: 20)
# Output: 
#   - se_cond_1_filtered: SummarizedExperiment with updated half-life segment identifiers
# ============================================


segment_HL_Rcpp_comp <- function(se_cond_1, se_cond_2, cores = 1, pen = 20) {
  # INCLUDE DECAY CONSTANT # PROBABLY MOVE TO OTHER PLACE IN CODE
  start_time <- Sys.time() 
  se_cond_1 <- decay_c_synth(se_cond_1)
  se_cond_2 <- decay_c_synth(se_cond_2)

  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  rowRanges(se_cond_1)$half_life_dif <- rowRanges(se_cond_2)$half_life-rowRanges(se_cond_1)$half_life
  rowRanges(se_cond_1)$log2_half_life_fold_change <- log2(rowRanges(se_cond_2)$half_life/rowRanges(se_cond_1)$half_life)
  rowRanges(se_cond_1)$decay_constant_1 <- rowRanges(se_cond_1)$decay_constant
  rowRanges(se_cond_1)$decay_constant_2 <- rowRanges(se_cond_2)$decay_constant
  rowRanges(se_cond_1)$synthesis_rate_1 <- rowRanges(se_cond_1)$synthesis_rate
  rowRanges(se_cond_1)$synthesis_rate_2 <- rowRanges(se_cond_2)$synthesis_rate
  rowRanges(se_cond_1)$intensity_2 <- rowRanges(se_cond_2)$intensity
  # the dataframe is sorted by strand and position.
  se_cond_1 <- inp_order(se_cond_1)
  #make the tmp_df using delay segments not position
  tmp_df <- inp_df(se_cond_1, "ID", "half_life_dif", "position_segment")
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  # the segmentation is performed on the position_segments, independent on
  # if they are (terminal) outliers or NAs.
  tmp_df[, "position_segment"] <- gsub("_O|_NA", "", tmp_df$position_segment)
  
  # although _NA will most probably will be dismissed here, because there
  # should be no half-life, if there is no delay
  tmp_df <- na.omit(tmp_df)
  
  unique_seg <- unlist(unique(tmp_df$position_segment))
  
  count <- 1
  
  pen <- pen
  # II. Dynamic Programming: the scoring function is interpreted

  segs <- mclapply(unique(tmp_df[["position_segment"]]), function(seg) {
    tmp_df[tmp_df[["position_segment"]] == seg, ]
  }, mc.cores = cores)
  
  Rcpp_seg_df <- data.frame(ID = numeric(0), 
                            half_life_dif = numeric(0), 
                            position_segment = character(0), 
                            FLT = numeric(0), 
                            strand = character(0), 
                            position = numeric(0))
  number = 0
  Rcpp_process_HL <- function(segment, pen){
    # Temporary dataframe containing the data of the current segment
    tmp_df <- data.frame(segment)
    # Temporary dataframe containing the Rcpp segmented data
    Rcpp_tmp_df <- solve_for_partition(1:nrow(tmp_df), tmp_df$half_life_dif, penalty = pen,min_seg = 3)
    # Extract the segment limits
    seg_list <- Rcpp_tmp_df$x
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
      Rcpp_tmp_df_subset$HL_dif_segment <- paste0("Dc_", number)
      # Return the subset data frame for the current segment
      Rcpp_tmp_df_subset
    })
    do.call(rbind, result)
  }
  
  final_segs <- lapply(segs, Rcpp_process_HL, pen)
  final_segs_df <- do.call(rbind, final_segs)
  matching_rows <- rownames(se_cond_1) %in% rownames(final_segs_df)
  # Subset the inp object based on matching rows
  se_cond_1_filtered <- se_cond_1[matching_rows, ]
  mcols(rowRanges(se_cond_1_filtered))[["HL_dif_segment"]] <- final_segs_df$HL_dif_segment
  rowRanges(se_cond_1_filtered)$HL_dif_segment <- final_segs_df$HL_dif_segment
  end_time <- Sys.time() 
  elapsed_t <- end_time - start_time
  print(elapsed_t)
  se_cond_1_filtered
}


