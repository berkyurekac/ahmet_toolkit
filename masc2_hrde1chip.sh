#!/bin/bash

#SBATCH -J macs2_pme1_inhrde1
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-48:00
#SBATCH -p IACT
#SBATCH --mail-user=ab2329@cam.ac.uk


macs2 callpeak -t /mnt/home1/miska/ab2329/chip/aligned/HRDE1^Q5347_aa188^E^N2^YA_raw^NA^NA_AA993^Cac77683.txt.sorted.bam  /mnt/home1/miska/ab2329/chip/aligned/HRDE1^Q5347_aa191^E^N2^YA_raw^NA^NA_AA994^C8977678.txt.sorted.bam -c /mnt/home1/miska/ab2329/chip/aligned/ce11_merged_input.bam -B --nomodel --extsize 500 --pvalue 0.001 --SPMR -g ce -n HRDE1^Q5347_aa191^E^N2^YA_peaks