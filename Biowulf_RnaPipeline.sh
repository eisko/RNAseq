#!/bin/bash

#To run on sbatch type:
# sbatch Biowulf_RnaPipeline.sh --cpus-per-task=8 --mem=20g 
#SBATCH --job-name EI_01_190102
 
#Load needed module
module load fastqc
module load hisat
module load samtools
module load stringtie
module load trimmomatic





#OUTPUT
#.txt file containing info on code run
#output files (fastq, sam, bam files...)
start_script=$(date +%s)


#DEFINED PATHS AND VARIABLES TO RUN CODE/SAMPLE

#Path to 'output' directory -> in data directory?
RnaSeq_home=/data/$USER/Vanc_RnaSeq


#merge files if needed - recommend doing this once in terminal before start
# Pipeline in case need to rerun part or whole script
# cat ~/code/RnaSeq/data/MUT3_5/MUT3_5_USR18003044L_HFC32DMXX_L2_1.fq.gz >> ~/code/RnaSeq/data/MUT3_5/MUT3_5_USR18003044L_HFC32DMXX_L1_1.fq.gz
# cat ~/code/RnaSeq/data/MUT3_5/MUT3_5_USR18003044L_HFC32DMXX_L2_2.fq.gz >> ~/code/RnaSeq/data/MUT3_5/MUT3_5_USR18003044L_HFC32DMXX_L1_2.fq.gz

#Path to file
file1=/data/$USER/Vanc_data/EI_01_1.fq.gz
file2=/data/$USER/Vanc_data/EI_01_2.fq.gz

#goto file name
#note: Af = affected, Uf = unaffected; Un = untreated, Vc = vancomycin; F = female, M = male
filename=EI_01


#Path to refs genome
#note: used GRCm38.p6 from NCBI
#see https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26
#used following command to download fasta file (.fna = Fasta of Nucleic Acids)
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRC38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#note: file is in .fna format (and can be renamed to .fa format?)
#May need to rename .fna to .fa
#ref=/data/$USER/ref/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz


#Path to hisat2 index (presumably downloaded)
#wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38.tar.gz
#note: hisat intakes genome index by Path/index_basename (leave off .#.ht extension!)
# use from Biowulf resources
hisat_index=$HISAT_INDEXES/grcm38/genome

#Path to gtf file
#used the following to download gtf/gff
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz
# note: stringtie doesn't accept gff.gz format as dowloaded from NCBI
# Use gunzip /PATH/genome.gff.gz to unzip file (gff file not too big)
gtf=/data/$USER/ref/Mus_musculus.GRCm38.94.gtf


#set-up output file
cd $RnaSeq_home
touch ${RnaSeq_home}/${filename}_pipeline_info.txt
out_file=${RnaSeq_home}/${filename}_pipeline_info.txt
echo "Pipeline was run for sample ${filename}">>$out_file
echo "Pipeline was run on $(date)">>$out_file
echo "Runtimes (in minutes) for pipeline:">>$out_file
#note: for each step of pipeline, start and end times (in seconds) are recorded and run time for that program
#and is appended to out_file

echo "starting fastqc"
start_fastqc=$(date +%s)
#Run fastQ? - already done with novogene
mkdir -p ${RnaSeq_home}/FastQC
fastqc $file1
fastqc $file2
#note: fastqc outputs files in same directory where input file located
mv ${RnaSeq_home}/Vanc_data/${filename}/*qc* ${RnaSeq_home}/FastQC

end_fastqc=$(date +%s)
runtime_fastqc=$((${end_fastqc}-${start_fastqc}))
echo "FastQC: $(bc -l <<< "${runtime_fastqc}/60" | cut -c 1-5)">>$out_file
echo "End of Fastqc"


# #trimmomatic (get rid of unpaired reads and adaptor sequences)
# mkdir -p ${RnaSeq_home}/trimmed
# cd ${RnaSeq_home}/trimmed
# trimmomatic PE -summary ${filename}_trimmomatic_sum.txt \
# $file1 \ #input file1
# $file2 \ #input file2
# ..../trimmed/${filename}_1P \ #output paired reads
# ../trimmed/${filename}_1U \ #output unpaired reads
# ../trimmed/${filename}_2P \
# ../trimmed/${filename}_2U \
# #Need to remove adaptors or novogene do this? (check with fastqc)
# ILLUMINACLIP:../../refs/illumina_multiplex.fa:2:30:10 \
# LEADING:3 \
# TRAILING:3 \
# SLIDINGWINDOW:4:20 \
# MINLEN:70

#ALIGNMENT with Hisat2

#build hisat2 index (optional- can download index from online)
#hisat2-build (ref_genome.fa) ~/PATH/genome_index

#note: Griffith tutorial uses other arguments to specify read groups etc., not sure why
#needed or why helpful

echo "Starting hisat2"

start_hisat2=$(date +%s)
#align with hisat2
#note- -dta flag needed for downstream assembly (ie for stringtie input)
mkdir -p ${RnaSeq_home}/hisat2
cd ${RnaSeq_home}/hisat2
hisat2 --rna-strandness RF --dta -p $SLURM_CPUS_PER_TASK  -x $hisat_index \
-1 $file1 -2 $file2 \
-S /lscratch/$SLURM_JOB_ID/${filename}.sam \
--summary-file ${filename}_alignStats.txt

end_hisat2=$(date +%s)
runtime_hisat2=$((${end_hisat2}-${start_hisat2}))
echo "hisat2: $(bc -l <<< "${runtime_hisat2}/60" | cut -c 1-5)">>$out_file

echo "end of hisat2"

echo "Start of samtools"
start_samtools=$(date +%s)
#Sort alignments by chromosomal coordinates and convert to bam format
mkdir -p ${RnaSeq_home}/sorted
cd ${RnaSeq_home}/sorted
samtools sort -@ $SLURM_CPUS_PER_TASK -o ${filename}_sorted.bam ${RnaSeq_home}/hisat2/${filename}.sam
end_samtools=$(date +%s)
runtime_samtools=$((${end_samtools}-${start_samtools}))
echo "Samtools: $(bc -l <<< "${runtime_samtools}/60" | cut -c 1-5)">>$out_file


# delete .sam files from hisat to free up space
rm /lscratch/$SLURM_JOB_ID/${filename}.sam

echo "end of samtools"

# #(optional) merge files of same group using picard (preserves heading titles?)
# #NOTE: only work if process multiple files at once
# #picard OUTPUT=outputfile.bam INPUT=inputfile_rep1.bam INPUT=inputfile_rep2.bam INPUT=inputfile_rep3.bam


echo "start of stringtie"
start_stringtie=$(date +%s)
#ASSEMBLY WITH stringtie
mkdir -p ${RnaSeq_home}/stringtie/$filename
cd ${RnaSeq_home}/stringtie/$filename
stringtie -p $SLURM_CPUS_PER_TASK -G $gtf -e -B -o ${filename}_transcripts.gtf -A ${filename}_gene_abundances.tsv ${RnaSeq_home}/sorted/${filename}_sorted.bam

end_stringtie=$(date +%s)
runtime_stringtie=$((${end_stringtie}-${start_stringtie}))
echo "Stringtie: $(bc -l <<< "${runtime_stringtie}/60" | cut -c 1-5)">>$out_file
echo "end of stringtie"



echo "Start of samtools"
start_samtools=$(date +%s)
#Index bam files for visualization (e.g. for IGV)
cd ${RnaSeq_home}/sorted
samtools index -@ $SLURM_CPUS_PER_TASK ${filename}.bam 
end_samtools=$(date +%s)
runtime_samtools=$((${end_samtools}-${start_samtools}))
echo "Samtools Index: $(bc -l <<< "${runtime_samtools}/60" | cut -c 1-5)">>$out_file


end_script=$(date +%s)
runtime_script=$((${end_script}-${start_script}))
echo "Script: $(bc -l <<< "${runtime_script}/60" | cut -c 1-5)">>$out_file

echo "end of script"
