# python script for generating swarm files
# needed to run file:
#   - PROJECT: path to folder containing all folders (hisat, sorted, stringtie...) needed to save files, see swarm outline for description of proper folder file setup 
#   - SAMPLES: samples.csv file containing sampleID to label subsequent created files,
#       and paths to fastq files (forward and reverse, i.e. there should be 2 fastq files per sample)
#   - need to provide path to reference files
#       - INDEX: hisat2 index files
#       - ANNOTATION: gff gene annotation files for stringtie

import pandas as pd
import numpy as np

# path to project folder containing other needed folders
PROJECT = "/data/$USER/RnaSeq"

# create pandas dataframe with phenotypic data
SAMPLES = pd.read_csv('samples.csv', sep=',')

# can download hisat index libraries from hisat2 site or download from biowulf reference folder
# may need to move to lscratch if running many jobs/processing many files/samples?
# can also build own indeces by running hisat2-build functions
INDEX = "/data/$USER/ref/grcm38/genome"
# include path to gtf file used to assign gene names to reads
ANNOTATION = "/data/$USER/ref/Mus_musculus.GRCm38.94.gtf"








# extract file names
sampleIDs = SAMPLES['SampleID']

# extract fastq file paths
fastqs1 = SAMPLES['Fastq_Path1']
fastqs2 = SAMPLES['Fastq_Path2']


############# SWARM 1 - FASTQC - run QC on Fastq files
# write first swarm script in pipeline - running QC on fastq files with fastqc
with open('fastq.swarm', 'w+') as file:
    # following to write as header in file
    file.write('# This is a swarm file for fastqc\n')
    file.write('# To run file:\n')
    file.write('# swarm -f PATH_TO_fastq.swarm --time 15:00 -p 2 -b 6 --module fastqc --job-name=fastqc_swarm\n')
    file.write('# NOTE: swarm script will make deposit fastqc files into same folder containing fastq files\n')
    file.write('# may need to move files in sinteractive if want fastqc files in a different folder')

    for path in fastqs1:
        file.write('fastqc %s' %path)
        file.write('\n')

    for path in fastqs2:
        file.write('fastqc %s' %path)
        file.write('\n')


###################### SWARM 2 - HISAT2 - Align sequences to reference genome
# write hisat2 swarm file
# NOTE: hisat2 command assumes starnded RF libraries (read1 = reverse, read2 = forward) and...
with open('hisat2.swarm', 'w+') as file:
    # following will write as header
    file.write('# To run:\n')
    file.write('# swarm -f PATH_TO_hisat2.swarm -g 5 -t 8 --module hisat --time 15:00 -b 4 --job-name=hisat2_swarm --merge-output\n')

    for sample, path1, path2 in zip(sampleIDs, fastqs1, fastqs2):
        # %(index, r1, r2, sampleID, sampleID)
        file.write('hisat2 --rna-strandness RF --dta -p $SLURM_CPUS_PER_TASK -x %s -1 %s -2 %s -S /scratch/$USER/%s.sam --summary-file %s/hisat2/%s_alignStats.txt\n' % (INDEX, path1, path2, sample, PROJECT, sample))


####################### SWARM 3 - SAMTOOLS - convert sam file to bam file
with open('samtools.swarm', 'w+') as file:
    # following will be file header
    file.write('# To run:\n')
    file.write('# swarm -f PATH_TO_samtools.swarm -t 4 -g 5 --time 15:00 -b 4 --job-name=samtools_swarm --merge-output --module samtools --dependency=afterok:hisat_swarm_JOBID\n')

    for sample in sampleIDs:
        file.write('samtools sort -@ $SLURM_CPUS_PER_TASK -o %s/sorted/%s_sorted.bam /scratch/$USER/%s.sam\n' % (PROJECT, sample, sample))


###################### SWARM 4 - STRINGTIE - Assemble transcripts
with open('stringtie.swarm', 'w+') as file:
    # following will be file header
    file.write('# To run:\n')
    file.write('swarm -f PATH_TO_stringtie.swarm --time 25:00 -b 4 --job-name=stringtie_swarm --module stringtie --merge-output --dependency=afterok:samtools_swarm_JOBID\n')
    
    for sample in sampleIDs:
        file.write('cd %s/stringtie/%s; stringtie -p $SLURM_CPUS_PER_TASK -G %s -e -B -o %s_transcripts.gtf -A %s_gene_abundances.tsv %s/sorted/%s_sorted.bam\n' % (PROJECT, sample, ANNOTATION, sample, sample, PROJECT, sample))



####################### SWARM 5 - SAMTOOLS INDEX - make bai files for viewing in IGV or other genome viewer applications
with open('samtools_index.swarm', 'w+') as file:
    # following will be header
    file.write('# To run:\n')
    file.write('# swarm -f ../myswarms/samtools_index.swarm --time=05:00 -b 4 --merge-output --job-name=samtools_index --module samtools --dependency=afterok:stringtie_swarm_JOBID\n')

    for sample in sampleIDs:
        file.write('cd %s/sorted; samtools index -@ $SLURM_CPUS_PER_TASK %s_sorted.bam\n' % (PROJECT, sample))
