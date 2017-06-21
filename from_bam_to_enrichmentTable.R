source("create_enrichment_table.R")
source("get_region_read_count_ratio.R")

args <- commandArgs(trailingOnly = TRUE)
sam <- args[1]
con <- args[2]
outname <- args[3]
gff <- "../Genomic_Features/saccharomyces_cerevisiae_R64-2-1_20150113.gff"

ratio <- get_region_read_count_ratio(sam, con, gff)

sambam <- open(BamFile(sam))
conbam <- open(BamFile(con))
create_enrichment_table(sambam, conbam, outname)