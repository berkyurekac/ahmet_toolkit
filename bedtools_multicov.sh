#!/bin/bash

#SBATCH -J bedtools_multicov
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-48:00
#SBATCH -p IACT
#SBATCH -o /mnt/home1/miska/ab2329/RNA-seq/pichip/bedtools_name.out
#SBATCH -e /mnt/home1/miska/ab2329/RNA-seq/pichip/bedtools_name.err
#SBATCH --mail-user=ab2329@cam.ac.uk


for file in *.fastq.gz_totalRNA_uniqueAligned.sortedByCoord.out.bam
do
NAME=$(basename "$file" .fastq.gz_totalRNA_uniqueAligned.sortedByCoord.out.bam)
bedtools coverage -S -a ce11_all_genes.bed -b "$file" -bed -counts > "$NAME".allgenes_expression_set2532unique.txt
done