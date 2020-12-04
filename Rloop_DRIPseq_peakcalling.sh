#!/bin/bash

#SBATCH -J macs2
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-48:00
#SBATCH -p IACT
#SBATCH --mail-user=ab2329@cam.ac.uk




macs2 callpeak -t SRR3993995Aligned.sortedByCoord.out.bam -c SRR3993993Aligned.sortedByCoord.out.bam -B --nomodel --extsize 500 --SPMR -g ce -n Rloop_wt