# ============================================
# prepare_microarray_fit_data
# Input:
#   - normalized_expression_iron.txt: Tab-delimited file containing iron condition expression data with columns:
#       - `start`, `end`: Genomic coordinates
#       - `orientation`: Strand orientation (+/-)
#       - `DFB*` columns: Expression values for iron condition replicates
#   - normalized_expression_standard.txt: Tab-delimited file containing standard condition expression data with columns:
#       - `start`, `end`: Genomic coordinates
#       - `orientation`: Strand orientation (+/-)
#       - `WT*` columns: Expression values for standard condition replicates
# Output:
#   - se_norm_iron: A SummarizedExperiment object containing iron condition expression matrix, metadata (colData), and genomic ranges (rowRanges)
#   - se_norm_standard: A SummarizedExperiment object containing standard condition expression matrix, metadata (colData), and genomic ranges (rowRanges)
# ============================================
# Laden der txt Dateien 
library(data.table)
library(readxl)
library(readr)
norm_iron_dt <- as.data.table(read_delim("~/Desktop/normalized_microarray_data_for_fits/normalized_expression_iron.txt"))
norm_iron_dt <- norm_iron_dt %>%
  filter(!is.na(start)) %>%
  mutate(position_1 = start,
         position_2 = end,
         strand = ifelse(orientation == "+", "-", "+"))

assay_data <- as.matrix(norm_iron_dt[, grep("^DFB", names(norm_iron_dt)), with = FALSE])
rownames(assay_data) <- norm_iron_dt[[1]]
time_points <- rep(c(0, 2, 4, 8, 16, 32, 64), times = 2)  # Two replicates
replicates <- rep(1:2, each = 7)
# Create colData
col_data <- DataFrame(
  timepoint = time_points,
  replicate = replicates
)
colnames(assay_data) <- paste0("T", time_points, "_R", replicates)
rowRanges <- GRanges(
  seqnames = Rle(rep("chr", nrow(norm_iron_dt))), 
  ranges = IRanges(
    start = as.numeric(norm_iron_dt$position_1),  # Convert start positions to numeric
    end = as.numeric(norm_iron_dt$position_2)    # Convert end positions to numeric
  ),
  strand = Rle(norm_iron_dt$strand),  # Strand information
  position = as.integer(Rle(norm_iron_dt$position_2))
)
se_norm_iron <- SummarizedExperiment(
  assays = list(counts = assay_data),
  colData = col_data,
  rowRanges = rowRanges
)



norm_standard_dt <- as.data.table(read_delim("~/Desktop/normalized_microarray_data_for_fits/normalized_expression_standard.txt"))

norm_standard_dt <- norm_standard_dt %>%
  filter(!is.na(start)) %>%
  mutate(position_1 = start,
         position_2 = end,
         strand = ifelse(orientation == "+", "-", "+"))

assay_data <- as.matrix(norm_standard_dt[, grep("^WT", names(norm_standard_dt)), with = FALSE])
rownames(assay_data) <- norm_standard_dt[[1]]
time_points <- c(0, 2, 4, 8, 16, 32, 64, 0, 2, 4, 8, 16, 32, 64, 0, 8, 2, 4)  # Two/Three replicates
replicates <- c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3)
# Create colData
col_data <- DataFrame(
  timepoint = time_points,
  replicate = replicates
)
colnames(assay_data) <- paste0("T", time_points, "_R", replicates)
rowRanges <- GRanges(
  seqnames = Rle(rep("chr", nrow(norm_standard_dt))), 
  ranges = IRanges(
    start = as.numeric(norm_standard_dt$position_1),  # Convert start positions to numeric
    end = as.numeric(norm_standard_dt$position_2)    # Convert end positions to numeric
  ),
  strand = Rle(norm_standard_dt$strand),  # Strand information
  position = as.integer(Rle(norm_standard_dt$position_2))
)
se_norm_standard <- SummarizedExperiment(
  assays = list(counts = assay_data),
  colData = col_data,
  rowRanges = rowRanges
)



