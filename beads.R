library(rbeads)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(IRanges)
library(digest)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)
library(methods)
library(devtools)
library(TxDb.Celegans.UCSC.ce11.refGene)

setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11")
# Get the paths of bam files
SX1316_1 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/beads_normalisation/mapped/SX1316_ChIP_rep1.sorted.bam'
a <- readGAlignments(SX1316_1)
b <- coverage(a)
c <- as.vector(b[8])

SX1316_2 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/Sx1316_ChIP_rep2.sorted.bam'

SX1316_3 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/SX1316_ChIP_ce11_rep3.sorted.bam'

SX1984_1 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11_rep2.sorted.bam'

SX1984_2 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11.sorted.bam'

SX1984_3 <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/mj261/SX1984_ChIP_ce11_rep3.sorted.bam'

# prepare merged input in bw format

input_bam <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/beads_normalisation/mapped/SX1316_input_rep1.sorted.bam'
d <- readGAlignments(input_bam)

map_bw <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/beads_normalisation/ce11.bw'
ref_fa <- '/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/beads_normalisation/Caenorhabditis_elegans.WBcel235ChIP_new.fa'


mappability <- import(map_bw)
mappability2 <- import(map_bw2)
mappability3 <- import(input_bam)

info <- seqinfo(mappability)
print(info)

info2 <- seqinfo(mappability2)
print(info2)

info3 <- seqinfo(mappability3)
print(info3)

info_bam <- seqinfo(SX1316_1)
# Set the directory where the output files will be crated
setwd(tempdir())

# Run BEADS for BAM input file
beads(SX1316_1, input_bam, map_bw, ref_fa)

# Run BEADS for SummedInput (BigWig) input file
beads(SX1316_1, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 200L, mapq_cutoff = 10L)

beads(SX1316_2, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 50L, mapq_cutoff = 10L, export = "SX1316_rep1")

beads(SX1316_2, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 50L, mapq_cutoff = 10L, export = "SX1316_rep1")

beads(SX1984_1, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 50L, mapq_cutoff = 10L, export = "SX1316_rep1")

beads(SX1984_2, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 50L, mapq_cutoff = 10L, export = "SX1316_rep1")

beads(SX1984_3, input_bam, map_bw, ref_fa, uniq = FALSE, insert = 50L, mapq_cutoff = 10L, export = "SX1316_rep1")

## Not run: 
## Run BEADS for BSgenome package, the reference genome package have to be installed prior to running this example
# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Celegans.UCSC.ce10")
# library(BSgenome.Celegans.UCSC.ce10)
beads(sample_bam, SummedInput_bw, map_bw, genome='ce10')

## Run BEADS for all BAM files in the directory
#lapply(dir(pattern='bam$'), beads, control=input, mappability=map_bw, genome=ref_fa)

## End(Not run)

## Not run: 
## In parallel:
# library(parallel)
# mclapply(dir(pattern='.bam$'), beads, control=input, mappability=map_bw, genome='ce10', mc.cores=parallel::detectCores())

## End(Not run)