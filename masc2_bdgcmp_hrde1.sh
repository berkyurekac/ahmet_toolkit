#!/bin/bash

#SBATCH -J macs2
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-48:00
#SBATCH -p IACT
#SBATCH --mail-user=ab2329@cam.ac.uk


macs2 bdgcmp -t SRR1054262_peaks_NOSPMR_treat_pileup.bdg  -c SRR1054262_peaks_NOSPMR_control_lambda.bdg -o SRR1054262_linear.bdg -m FE

macs2 bdgcmp -t SRR1054263_peaks_NOSPMR_treat_pileup.bdg  -c SRR1054263_peaks_NOSPMR_control_lambda.bdg -o SRR1054263_linear.bdg -m FE

macs2 bdgcmp -t SRR1054264_peaks_NOSPMR_treat_pileup.bdg  -c SRR1054264_peaks_NOSPMR_control_lambda.bdg -o SRR1054264_linear.bdg -m FE

