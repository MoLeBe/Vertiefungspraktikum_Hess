# ============================================
# Required Libraries
# ============================================
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
library("gridExtra")
library("fuzzyjoin")
library("data.table")


# ============================================
# detect_outliers_cleaned
# Input: numeric vector of column data
# Output: list with:
#   - cleaned_data: outliers replaced with NA
#   - outliers: logical vector of outlier positions
# ============================================
detect_outliers_cleaned <- function(column_data) {
  median_value <- median(column_data, na.rm = TRUE)
  mad_value <- mad(column_data, constant = 1.4826, na.rm = TRUE)
  
  if (is.na(mad_value) || mad_value == 0) mad_value <- 1
  
  modified_z_score <- 0.6745 * (column_data - median_value) / mad_value
  outliers <- abs(modified_z_score) > 3.5
  cleaned_data <- ifelse(outliers, NA, column_data)
  
  list(cleaned_data = cleaned_data, outliers = outliers)
}

# ============================================
# decay_c_synth
# Input: SummarizedExperiment with 'half_life' and 'intensity' in rowRanges
# Output: Same SummarizedExperiment with added 'decay_constant' and 'synthesis_rate'
# ============================================
decay_c_synth <- function(inp){
  rowRanges(inp)$decay_constant <- log(2)/rowRanges(inp)$half_life
  rowRanges(inp)$synthesis_rate <- rowRanges(inp)$intensity * rowRanges(inp)$decay_constant
  return(inp)
}

# ============================================
# preprocess_comp_se
# Input:
#   - comp_se: SummarizedExperiment with differential segmentation data
#   - strand_type: "+" / "-" / "both"
# Output: data.frame of filtered, minimal representation
# ============================================
preprocess_comp_se <- function(comp_se, strand_type) {
  rowRanges(comp_se) %>%
    as.data.frame() %>%
    filter(if (strand_type == "both") TRUE else strand == strand_type) %>%
    mutate(gene = NA) %>%
    select(
      half_life_dif, log2_half_life_fold_change, HL_dif_segment,
      log2_intensity_fold_change, intensity_log2_fc_segment, position,
      gene, decay_constant_1, decay_constant_2, synthesis_rate_1,
      synthesis_rate_2, ID, intensity, intensity_2, strand, ID
    )
}

# ============================================
# match_genes
# Input:
#   - df: dataframe with positions and strand info
#   - gff3_data: gene annotation (start, end, gene_locus, strand)
# Output: df with gene loci matched based on position and strand
# ============================================
match_genes <- function(df, gff3_data) {
  df_dt <- as.data.table(df)
  gff3_dt <- as.data.table(gff3_data)
  setkey(gff3_dt, start, end)
  
  matched <- foverlaps(
    x = df_dt[, .(query_start = position, query_end = position, ID, strand)],
    y = gff3_dt[, .(start, end, gene_locus, strand)],
    by.x = c("query_start", "query_end"),
    by.y = c("start", "end"),
    type = "within",
    nomatch = 0
  )
  
  matched <- matched[strand == i.strand]
  df_dt <- merge(df_dt, matched[, .(query_start, gene_locus)], 
                 by.x = "position", by.y = "query_start", 
                 all.x = TRUE)
  df_dt[, gene := gene_locus]
  return(as.data.frame(df_dt))
}

# ============================================
# group_by_seg
# Input: dataframe with per-position expression data
# Output: summary dataframe grouped by segment with cleaned + aggregated values
# ============================================
group_by_seg <- function(df) {
  df %>%
    group_by(intensity_log2_fc_segment) %>%
    mutate(
      across(
        c(half_life_dif, log2_half_life_fold_change,
          log2_intensity_fold_change, decay_constant_1,
          decay_constant_2, synthesis_rate_1, synthesis_rate_2,
          intensity, intensity_2),
        ~ detect_outliers_cleaned(.x)$cleaned_data,
        .names = "{.col}_clean"
      )
    ) %>%
    summarise(
      across(ends_with("_clean"), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
      gene = paste(setdiff(unique(gene), "NULL"), collapse = ", "),
      position_1 = min(position, na.rm = TRUE),
      position_2 = max(position, na.rm = TRUE),
      mean_log_decay_constant_fold_change = log(mean(decay_constant_1_clean, na.rm = TRUE) / mean(decay_constant_2_clean, na.rm = TRUE)),
      mean_log_synthesis_rate_fold_change = log(mean(synthesis_rate_2_clean, na.rm = TRUE) / mean(synthesis_rate_1_clean, na.rm = TRUE)),
      mean_intensity = mean(intensity, na.rm = TRUE),
      mean_intensity_2 = mean(intensity_2, na.rm = TRUE),
      IDs = paste(unique(ID), collapse = ", ")
    ) %>%
    ungroup()
}

# ============================================
# get_dif_se
# Input:
#   - se_cond_1: SummarizedExperiment for condition 1
#   - se_cond_2: SummarizedExperiment for condition 2
#   - gff3_file_path: path to GFF3 annotation file
# Output: list of sense/antisense/full expression summaries and top changes
# ============================================
get_dif_se <- function(se_cond_1, se_cond_2, gff3_file_path){
  start_time <- Sys.time()
  
  comp_se <- se_cond_1 %>%
    segment_HL_Rcpp_comp(se_cond_2) %>%
    segment_int_Rcpp_comp(se_cond_2)
  
  sense_tmp_comp_df <- preprocess_comp_se(comp_se, "+")
  antisense_tmp_comp_df <- preprocess_comp_se(comp_se, "-")
  full_tmp_df <- preprocess_comp_se(comp_se, "both")
  
  gff3_data <- gff3_preprocess(gff3_file_path)
  gff3_sense <- gff3_data %>% filter(strand == "+")
  gff3_antisense <- gff3_data %>% filter(strand == "-")
  
  sense_tmp_comp_df <- match_genes(sense_tmp_comp_df, gff3_sense)
  antisense_tmp_comp_df <- match_genes(antisense_tmp_comp_df, gff3_antisense)
  full_tmp_comp <- match_genes(full_tmp_df, gff3_data)
  
  grouped_sense_tmp_comp_df <- group_by_seg(sense_tmp_comp_df)
  grouped_antisense_tmp_comp_df <- group_by_seg(antisense_tmp_comp_df)
  grouped_full_comp_df <- bind_rows(sense_tmp_comp_df, antisense_tmp_comp_df) %>% group_by_seg()
  
  top_sense_comp_df <- great_dif(grouped_sense_tmp_comp_df)
  top_antisense_comp_df <- great_dif(grouped_antisense_tmp_comp_df)
  
  end_time <- Sys.time()
  elapsed_t <- end_time - start_time
  print(elapsed_t)
  
  return(list(
    sense_differential_expression = grouped_sense_tmp_comp_df,
    antisense_differential_expression = grouped_antisense_tmp_comp_df,
    full_differential_expression = grouped_full_comp_df,
    top_sense_differential_expression = top_sense_comp_df,
    top_antisense_differential_expression = top_antisense_comp_df
  ))
}
