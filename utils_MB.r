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
library("showtext")

char_TI <- function(y, hig = 0.075) {
  y <- y / y[1] # normalizes the input
  tmp <-
    max(y, na.rm = TRUE) - y[1] # subtracts t0 point from the highest point
  res <- -1 # not TI exept..
  if (tmp > hig) {
    # ..the difference is higher than 7.5%
    res <- 1 # if yes, the output is 1
  }
  res
}

relativity <- function(x) {
  (x / max(x))
}

finding_above_backg <- function(df, bg) {
  res <- rep("_ABG_",nrow(df))
  for (j in seq_len(nrow(df))) {
    tmp_last <- df[j, ncol(df)]
    if (!is.na(tmp_last)) {
      if (tmp_last < bg) {
        res[j] <- "_"
      }
    }
  }
  res
}

assert <- function(expr, error) {
  if (!expr) {
    stop(error, call. = FALSE)
  }
}

encode_FLT <- function(obj, rows, rep){
  bi <- intToBits(0)
  bi[rep] <- as.raw(1)
  int<-packBits(bi,"integer")
  rowRanges(obj)$FLT[rows] <- bitwOr(rowRanges(obj)$FLT[rows], int)
  obj
}


decode_FLT <- function(obj){
  #get int vector
  int_vec <- rowRanges(obj)$FLT
  #get logical from bits
  logi <- lapply(int_vec, function(x){as.logical(intToBits(x))})
  #representation fo replicates to filter
  repli <- lapply(logi, function(x){metadata(obj)$replicates[x]})
  #logical of correct length
  logi2 <- lapply(repli, function(x){colData(obj)$replicate %in% x})
  #logical matrix
  log_mat <- do.call(rbind,logi2)
  log_mat
}

inp_order <- function(inp){
  #the columns are ordered in a way, that t0 is in first position
  #the rows are order by strand and position increasing(plus) decreasing(minus)
  assert(all(c("timepoint", "replicate") %in% names(colData(inp))),
         "timepoint and replicate must be columns in the colData")
  assert("position" %in% names(mcols(rowRanges(inp))),
         "position must be columns in the rowRanges")
  ord_col <- order(colData(inp)$replicate,colData(inp)$timepoint)
  ord_row <- order(strand(inp))
  inp <- inp[ord_row, ord_col]
  ord_plus <- order(rowRanges(inp[strand(inp) == "+",])$position)
  ord_minus <- order(rowRanges(inp[strand(inp) == "-",])$position,
                     decreasing = TRUE)
  inp[strand(inp) == "+",] <- inp[strand(inp) == "+",][ord_plus,]
  inp[strand(inp) == "-",] <- inp[strand(inp) == "-",][ord_minus,]
  inp
}

inp_df <- function(inp, ...){
  cols <- c(..., "FLT")
  assert(all(cols %in% names(mcols(rowRanges(inp)))),
         paste0("one of the necessary columns in rowRanges is missing: ",
                paste0(cols,collapse = ", ")))
  #get the integer representing all replicates
  bi<-intToBits(0)
  bi[metadata(inp)$replicate] <- as.raw(1)
  int<-packBits(bi,"integer")
  #get the data frame
  tmp_df <- rowRanges(inp)[,cols]
  tmp_df <- as.data.frame(mcols(tmp_df))
  tmp_df$strand <- decode(strand(inp))
  tmp_df$position <- rowRanges(inp)$position
  #filter the data frame
  tmp_df <- tmp_df[tmp_df$FLT != int,]
  tmp_df
}

tmp_df_rev <- function(tmp_df, stra){
  tmp_df[tmp_df$strand == stra, "position"] <-
    ((tmp_df[tmp_df$strand == stra,
             "position"]) - (tmp_df[tmp_df$strand == stra, ]
                             [nrow(tmp_df[tmp_df$strand == stra, ]),
                               "position"])) * -1
  tmp_df
}

round_5 <- function(x){
  ceiling(x/5)*5
}

round_10 <- function(x){
  ceiling(x / 10) * 10
}

assert <- function(expr, error) {
  # If the expression is FALSE, stop and display the error message
  if (!expr) {
    stop(error, call. = FALSE)
  }
}

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
        str_replace_all("^c\\(|\\)$", "") %>%   # Remove 'c(...)'
        str_replace_all("\\s+,", ",") %>%      # Remove spaces before commas
        str_replace_all("\\\\", "") %>%        # Remove backslashes
        str_trim()                             # Remove leading/trailing spaces
    )
  return(tmp)
}

theme_void_yax <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(family = "Oxanium", face = "bold"),
        text = element_text(family = "Oxanium"),
        legend.position="none",
        panel.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
}