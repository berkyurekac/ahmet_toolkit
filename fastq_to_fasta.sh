#!/bin/bash

#SBATCH -J fastq_to_fasta
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-48:00
#SBATCH -p 1804
#SBATCH -o /mnt/home1/miska/ab2329/RNA-seq/samtools_name.out
#SBATCH -e /mnt/home1/miska/ab2329/RNA-seq/samtools_name.err
#SBATCH --mail-user=ab2329@cam.ac.uk




for file in *_trimmed.fastq; do
    NAME=$(basename "$file" _trimmed.fastq)
    fastq_to_fasta -r -i "$file" -o "${NAME}"_trimmed.fa
done