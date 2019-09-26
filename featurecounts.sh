#!/bin/bash
# sbatch featurecounts.sh
#SBATCH --job-name feature_counts_trial

#Load needed module
module load subread

gtf=/data/$USER/Vanc_RnaSeq/ref/Mus_musculus.GRCm38.94.gtf
inputfile = /data/iskoec/Vanc_RnaSeq/sorted/EI_01_sorted.bam

cd /data/$USER/data/iskoec/Vanc_RnaSeq/featurecounts

# see: https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf
# for user userguide
# -a = annotation file (gtf)
# -p = paired end reads
# -T = number of threads to use
# -t = feature used to calculate/determine meta features
# -g = meat-feature to be counted
featureCounts -p -a $gtf -t exon -g gene_id -o EI01_counts.txt /data/iskoec/Vanc_RnaSeq/sorted/EI_01_sorted.bam
