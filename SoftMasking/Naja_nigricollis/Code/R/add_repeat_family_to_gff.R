library(rtracklayer)
library(dplyr)

# Read directory from command-line argument
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

landscape_file <- grep('landscape$', list.files(dir), value = TRUE)
gff_file <- grep('reformat.gff3$', list.files(dir), value = TRUE)

landscape <- read.table(paste0(dir, landscape_file), skip = 6, fill = TRUE)# remove first 6 rows
cutoff_row <- which(landscape$V1 == "Coverage")[1] # keep only the df before "Coverage" appears in the dataframe
landscape <- landscape[1:cutoff_row - 1, ]

gff <- rtracklayer::import(paste0(dir, gff_file))
gff <- as.data.frame(gff)

gff <- gff %>%
  separate(Target, into = c("temp", "Target"), sep = " ", extra = "merge", fill = "right", remove = FALSE) %>%
  left_join(landscape %>% select(V1, V2), by = c('temp' = 'V2'), relationship = "many-to-many") %>%
  relocate(V1, .before = "temp") %>%
  rename('Name' = V1) %>%
  mutate(Target = paste(temp, Target, sep = ' ')) %>%
  select(-temp)

rtracklayer::export(gff, con = paste0(dir, gsub('reformat', 'repeatFamily', gff_file)), format = "gff3")
