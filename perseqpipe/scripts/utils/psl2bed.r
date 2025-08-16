# Script to transform PSL file from BLAT in BED12
# By default, keep only EndToEnd alignments without mismatch from BLAT ->
# this is done by keeping only alignment where query length == alignment length

library(data.table)
library(dplyr)

args <- commandArgs(TRUE)

tab <- fread(args[1], skip = 5)

# keep just rows where query length (V11) is equal to blockSize (V19) == alignment
tab_fil <- tab %>%
  filter(V11 == gsub(",","",V19))

# create two columns that we need for BED12 (score and rgb)
tab_fil$V22 <- "0"
tab_fil$V23 <- "0,0,0"

# reorder filtered PSL file to create BED12
tab_bed <- tab_fil[,c("V14","V16","V17","V10","V22","V9","V16","V17","V23","V18","V19","V20")]

fwrite(tab_bed, args[2], sep="\t", col.names = F)
