# EBV seqlogo script (top50 -> drop '*' -> r3_r2_ratio >= 10 -> most common length = 48 peptides;
# identical to the EBV subset of Selected_Peptides.csv)
if (!require(ggseqlogo)) {
  install.packages("ggseqlogo")
  library(ggseqlogo)
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}

file_path <- "../data/top10000_combine_9mer_2311-4TCR-EBV-SCD-round123.csv"
if (!file.exists(file_path)) {
  stop(paste0("File not found: ", file_path))
}
df <- read.csv(file_path)

if (!"aa" %in% names(df)) {
  stop("The input table has no 'aa' column; check the column names.")
}

# top 50
peptide_sequences <- df %>% slice(1:50)

# drop '*' peptides
peptide_sequences <- peptide_sequences %>% filter(!str_detect(aa, "\\*"))
if (nrow(peptide_sequences) == 0) {
  stop("All sequences were filtered out (starred peptides).")
}

# keep r3_r2_ratio >= 10
if (!"r3_r2_ratio" %in% names(peptide_sequences)) {
  stop("The input table has no 'r3_r2_ratio' column; check the column names.")
}
peptide_sequences <- peptide_sequences %>% filter(r3_r2_ratio >= 10)
if (nrow(peptide_sequences) == 0) {
  stop("No sequences satisfy r3_r2_ratio >= 10.")
}

# extract the 'aa' column as a character vector
peptide_sequences <- peptide_sequences$aa %>% na.omit() %>% as.character()
if (length(peptide_sequences) == 0) {
  stop("No valid peptide sequences.")
}

# drop empty strings, if any
empty_seq_count <- sum(nzchar(peptide_sequences) == FALSE)
if (empty_seq_count > 0) {
  warning(paste0("Dropped ", empty_seq_count, " empty-string sequences."))
}
peptide_sequences <- peptide_sequences[nzchar(peptide_sequences)]

# keep the most common length
sequence_lengths <- nchar(peptide_sequences)
most_common_length <- names(sort(table(sequence_lengths), decreasing = TRUE)[1]) %>% as.numeric()
valid_sequences <- peptide_sequences[sequence_lengths == most_common_length]

p <- ggseqlogo(valid_sequences,
               method = "probability",
               stack_width = 1,
               font = "helvetica_bold",
               col_scheme = "chemistry") +
  theme_logo()

dir.create("../results/seqlogos", showWarnings = FALSE, recursive = TRUE)
ggsave("../results/seqlogos/bold_top200_9mer_EBV_mutation.pdf",
       plot = p, width = 8, height = 6, dpi = 300)
print(paste("EBV seqlogo done, peptides:", length(valid_sequences)))
