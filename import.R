library(rtracklayer)
library(GenomicRanges)
setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/AMA-1/ce11/wt/normalized_replicates")

bw_combine_reps <- function(r1=NULL, r2=NULL, outname=paste0('Combined_BigWig_', gsub('-|:| ', '', Sys.time()), '.bw')) {
  
  library(rtracklayer) 
  message('Using replicate with r1=', r1, ' and r2=', r2 )
  paths <- c(r1, r2)
  
  bwfl <- BigWigFileList(unlist(paths))
  names(bwfl) <- basename(paths)
  
  suppressWarnings(rm(list='bw'))
  for(n in names(bwfl)) {
    message('Summarizing ', n)
    if(exists('bw')) {
      bw <- bw + import.bw(bwfl[[n]], as='Rle')
    } else {
      bw <- import.bw(bwfl[[n]], as='Rle')
    }
  }
  bw <- bw / length(bwfl)
  
  
  message('Exporting: ', outname)
  export.bw(bw, outname)
  
  return( outname )
}

# combined log2
bw_combine_reps(r1="SX1984_rep1_Zscore.bw", r2="SX1984_rep3_Zscore.bw")

# combined log2
bw_combine_reps(r1="HRDE1^Q5347_aa188^E^N2^YA_BEADSNQNU^log2^ce11_AA993^B6277752.bw", r2="HRDE1^Q5347_aa191^E^N2^YA_BEADSQ10NU^log2^ce11_AA994^B7177745.bw")

bw_combine_reps(r1="SNPC4^AB290_ak01^E^YL457^YA_BEADSQ10NU^log2^ce11_AK009^Ba653010.bw", r2="SNPC4^AB291_ak02^E^YL457^YA_BEADSQ10NU^log2^ce11_AK010^B1a53009.bw")

bw_combine_reps(r1="TOFU4^ab290_akx1^E^N2^YA_BEADSQ10NU^log2^ce11_AK013^B3453012.bw", r2="TOFU4^NA_ak06^E^N2^YA_BEADSQ10NU^log2^ce11_AK041^B3457993.bw")

bw_combine_reps(r1="TOFU5^AB290_ak05^E^TOFU5-GFP^YA_BEADSQ10NU^log2^ce11_AK004^Bbb53001.bw", r2="TOFU5^NA_ak08^E^N2^YA_BEADSQ10NU^log2^ce11_AK045^B7557966.bw")

#combined zscores
bw_combine_reps(r1="PRDE1^ksoj_ak03^E^N2^YA_BEADSQ10NU^zscore^ce11_AK005^B4744894.bw", r2="PRDE1^ksoj_ak04^E^N2^YA_BEADSQ10NU^zscore^ce11_AK006^Bc244880.bw")

bw_combine_reps(r1="SNPC4^AB290_ak01^E^YL457^YA_BEADSQ10NU^zscore^ce11_AK009^B9d45030.bw", r2="SNPC4^AB291_ak02^E^YL457^YA_BEADSQ10NU^zscore^ce11_AK010^Bf445046.bw")

bw_combine_reps(r1="TOFU4^ab290_akx1^E^N2^YA_BEADSQ10NU^zscore^ce11_AK013^B3444899.bw", r2="TOFU4^NA_ak06^E^N2^YA_BEADSQ10NU^zscore^ce11_AK041^B5b57992.bw")

bw_combine_reps(r1="TOFU5^AB290_ak05^E^TOFU5-GFP^YA_BEADSQ10NU^zscore^ce11_AK004^Bda44816.bw", r2="TOFU5^NA_ak08^E^N2^YA_BEADSQ10NU^zscore^ce11_AK045^B5557965.bw")

#combined linear
bw_combine_reps(r1="PRDE1^ksoj_ak04^E^N2^YA_BEADSQ10NU^linear^ce11_AK006^B7e44878.bw", r2="PRDE1^ksoj_ak03^E^N2^YA_BEADSQ10NU^linear^ce11_AK005^B9344891.bw")

bw_combine_reps(r1="SNPC4^AB290_ak01^E^YL457^YA_BEADSQ10NU^linear^ce11_AK009^Bca45027.bw", r2="SNPC4^AB291_ak02^E^YL457^YA_BEADSQ10NU^linear^ce11_AK010^B4a45041.bw")

bw_combine_reps(r1="TOFU4^NA_ak06^E^N2^YA_BEADSQ10NU^linear^ce11_AK041^B2657991.bw", r2="TOFU4^ab290_akx1^E^N2^YA_BEADSQ10NU^linear^ce11_AK013^Bc844898.bw")

bw_combine_reps(r1="TOFU5^AB290_ak05^E^TOFU5-GFP^YA_BEADSQ10NU^linear^ce11_AK004^B3e44815.bw", r2="TOFU5^NA_ak08^E^N2^YA_BEADSQ10NU^linear^ce11_AK045^B6657964.bw")

