# Rebuilt script (2026-10-07): the original 1G4 seqlogo script was not archived;
# this script recreates it with the verified criteria.
# Criteria: top150 -> drop '*' -> r3_r2_ratio >= 4 -> most common length = 105 peptides
# (identical to the 1G4 subset of Selected_Peptides.csv; the output matches
#  results/seqlogos/bold_top200_9mer_1G4_mutation.pdf)
if (!require(ggseqlogo)) {
  install.packages("ggseqlogo")
  library(ggseqlogo)
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}

file_path <- "../data/top10000_combine_9mer_2311-4TCR-1G4-SCD-round123.csv"
if (!file.exists(file_path)) {
  stop(paste0("File not found: ", file_path))
}
df <- read.csv(file_path)

if (!"aa" %in% names(df)) {
  stop("The input table has no 'aa' column; check the column names.")
}

# top 150 -> drop '*' peptides -> r3_r2_ratio >= 4
peptide_sequences <- df %>% slice(1:150) %>%
  filter(!str_detect(aa, "\\*")) %>%
  filter(r3_r2_ratio >= 4)
if (nrow(peptide_sequences) == 0) {
  stop("No sequences left after filtering.")
}

peptide_sequences <- peptide_sequences$aa %>% na.omit() %>% as.character()
peptide_sequences <- peptide_sequences[nzchar(peptide_sequences)]

# most common length (105 peptides, no duplicates)
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
ggsave("../results/seqlogos/bold_top200_9mer_1G4_mutation.pdf",
       plot = p, width = 8, height = 6, dpi = 300)
print(paste("1G4 seqlogo done, peptides:", length(valid_sequences)))
