# TP53 Gene Analysis - Human GRCh38

# Install and load required packages
required_cran <- c("ggplot2", "dplyr")
required_bioc <- c("biomaRt", "GenomicRanges", "Gviz")

missing_cran <- required_cran[!required_cran %in% installed.packages()[, "Package"]]
if(length(missing_cran)) install.packages(missing_cran)

if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- required_bioc[!required_bioc %in% installed.packages()[, "Package"]]
if(length(missing_bioc)) BiocManager::install(missing_bioc)

library(biomaRt)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(Gviz)

# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# TP53 genomic region
chr <- "17"
start_pos <- 7661779
end_pos <- 7687546

# Gene info
gene_info <- getBM(
  attributes = c(
    'ensembl_gene_id', 'external_gene_name', 'chromosome_name',
    'start_position', 'end_position', 'strand', 'gene_biotype', 'description'
  ),
  filters = c('chromosome_name', 'start', 'end'),
  values = list(chr, start_pos, end_pos),
  mart = ensembl
)

print("Gene Info:")
print(gene_info)

# Transcript info
transcript_info <- getBM(
  attributes = c(
    'ensembl_transcript_id',
    'transcript_start',
    'transcript_end',
    'transcript_length',
    'transcript_biotype',
    'transcript_is_canonical'
  ),
  filters = 'external_gene_name',
  values = 'TP53',
  mart = ensembl
)

canonical <- transcript_info %>% filter(transcript_is_canonical == 1)
canonical_id <- canonical$ensembl_transcript_id[1]

print("Canonical Transcript:")
print(canonical)

# CDS info
cds_info <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "genomic_coding_start",
    "genomic_coding_end",
    "cds_start",
    "cds_end",
    "cds_length"
  ),
  filters = "ensembl_transcript_id",
  values = canonical_id,
  mart = ensembl
)

print("Canonical CDS Info:")
print(cds_info)

write.csv(cds_info, "TP53_canonical_CDS.csv", row.names = FALSE)

# Transcript visualization
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

ensemblTrack <- BiomartGeneRegionTrack(
  genome = "hg38",
  chromosome = chr,
  start = start_pos,
  end = end_pos,
  name = "TP53",
  biomart = ensembl,
  strand = "-"
)

plotTracks(
  list(itrack, gtrack, ensemblTrack),
  from = start_pos,
  to = end_pos,
  transcriptAnnotation = "symbol",
  background.title = "lightblue",
  col.title = "black",
  main = "TP53 Gene Structure (GRCh38)"
)

write.csv(transcript_info, "TP53_transcript_info.csv", row.names = FALSE)
