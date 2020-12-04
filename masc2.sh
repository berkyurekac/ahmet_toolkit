#!/bin/bash

#SBATCH -J macs2
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-48:00
#SBATCH -p IACT
#SBATCH --mail-user=ab2329@cam.ac.uk



for (( i = 62; i <= 64; i++ ))
  do
macs2 callpeak -t SRR10542${i}_1.sorted.bam -c SRR1054265_1.sorted.bam SRR1054266_1.sorted.bam -B --nomodel --extsize 500 -g ce -n SRR10542${i}_peaks_NOSPMR
done