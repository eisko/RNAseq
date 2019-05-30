#!:/bin/bash

# make sure to run script from /data/$USER to create project folders

# Name project:
project = 'Project_name'

# script to generate following folder structure:

mkdir /data/$USER/$project
cd /data/$USER/$project
# /data/$USER/PROJECT
# -Exp_data
#     |-contains raw fastq files
#     |-sample1_1.fq.gz
#     |-sample1_2.fq.gz
#     |-sample2_1.fq.gz
#     |
#     |
mkdir ref
# -ref
#     |-contains hisat indexes and annotation file (gtf)
#     |-grcm38 #folder w/ hisat indexes
#         |-genome.1.ht2
#         |-genome
#     |-Mus_musculus.GRCm38.94.gtf.gz
mkdir fastqc
# -fastqc
#     |-sample1_1.zip
#     |-sample1_1.html
#     |
#     |
mkdir hisat2
# -hisat2
#     |-sample1_alignStats.txt
#     |-sample2_alignStats.txt
mkdir sorted
# -sorted
#     |-sample1.bam
#     |-sample1.bam.bai
#     |-sample2.bam
#     |-sample2.bam.bai
mkdir stringtie # also need to create folders for every sample...
# -stringtie
#     |-sample1 #folder containing file output of stringtie, .ctab files for ballgown input
#         |-sample1_gene_abundances.tsv
#         |-sample1_transcripts.gtf
#         |-e_data.ctab
#         |-e2t.ctab
#         |-i_data.ctab
#         |-i2t.ctab
#         |-t_data.ctab
#     |-sample2

