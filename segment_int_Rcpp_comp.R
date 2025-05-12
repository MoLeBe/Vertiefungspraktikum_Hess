# ============================================
# segment_int_Rcpp_comp
# Input: 
#   - se_cond_1: SummarizedExperiment for condition 1
#   - se_cond_2: SummarizedExperiment for condition 2
#   - cores: number of cores for parallel processing (default: 1)
#   - pen: penalty value for dynamic programming (default: 20)
# Output: 
#   - se_cond_1_filtered: SummarizedExperiment with updated intensity segment identifiers
# ============================================

segment_int_Rcpp_comp <- function(se_cond_1,se_cond_2, cores = 1, pen = 20) {
  start_time <- Sys.time()
  # Ensure the data has the same rows
  common_rows <- intersect(rownames(se_cond_1), rownames(se_cond_2))
  se_cond_2 <- se_cond_2[common_rows, ]
  # I.Preparations: the dataframe is configured and some other variables are
  # assigned
  rowRanges(se_cond_1)$log2_intensity_fold_change <-  log2(rowRanges(se_cond_2)$intensity/rowRanges(se_cond_1)$intensity)
  
  # the dataframe is sorted by strand and position.
  se_cond_1 <- inp_order(se_cond_1)
  #make the tmp_df
  tmp_df <- inp_df(se_cond_1, "ID", "log2_intensity_fold_change", "HL_dif_segment")
  #revert the order in plus
  tmp_df <- tmp_df_rev(tmp_df, "-")
  
  # here it is important that (terminal) outliers and NAs are considered, as all
  # bins have an intensity that needs to be grouped.
  tmp_df[, "HL_dif_segment"] <- gsub("_O|_NA", "", tmp_df$HL_dif_segment)
  
  # this is just for safety, nothing should be omitted here
  tmp_df <- na.omit(tmp_df)
  
  
  unique_seg <- unlist(unique(tmp_df$HL_dif_segment))
  
  count <- 1
  pen <- pen
  # II. Dynamic Programming: the scoring function is interpreted
  
  
  
  segs <- mclapply(unique(tmp_df[["HL_dif_segment"]]), function(seg) {
    tmp_df[tmp_df[["HL_dif_segment"]] == seg, ]
  }, mc.cores = cores)
  
  
  Rcpp_seg_df <- data.frame(ID = numeric(0),
                            log2_intensity_fold_change = numeric(0),
                            HL_dif_segment = character(0),
                            FLT = numeric(0),
                            strand = character(0),
                            position = numeric(0))
  number = 0
  Rcpp_process_int <- function(segment, pen){
    # Temporary dataframe containing the data of the current segment
    tmp_df <- data.frame(segment)
    # Temporary dataframe containing the Rcpp segmented data
    Rcpp_tmp_df <- solve_for_partition(1:nrow(tmp_df), tmp_df$log2_intensity_fold_change, penalty = pen, min_seg = 3)
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
      Rcpp_tmp_df_subset$intensity_log2_fc_segment <- paste0("I_", number)
      # Return the subset data frame for the current segment
      Rcpp_tmp_df_subset
    })
    do.call(rbind, result)
  }
  
  final_segs <- lapply(segs, Rcpp_process_int, pen)
  final_segs_df <- do.call(rbind, final_segs)
  matching_rows <- rownames(se_cond_1) %in% rownames(final_segs_df)
  # Subset the inp object based on matching rows
  se_cond_1_filtered <- se_cond_1[matching_rows, ]
  mcols(rowRanges(se_cond_1_filtered))[["intensity_log2_fc_segment"]] <- final_segs_df$intensity_log2_fc_segment
  rowRanges(se_cond_1_filtered)$intensity_log2_fc_segment <- final_segs_df$intensity_log2_fc_segment
  end_time <- Sys.time() 
  elapsed_t <- end_time - start_time
  print(elapsed_t)
  se_cond_1_filtered
  


}
