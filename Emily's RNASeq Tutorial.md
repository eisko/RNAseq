# Emily's RNASeq Tutorial

This document is meant to be a brief overview of the tools and commands I used to analyze the samples. 
For a more a more comprehensive tutorial, I recommend the RNASeq tutorial maintained by the [**Griffith Lab**](https://rnabio.org/). 
For more information on individual tools, I recommend looking up the documentation/manuals from the developer's sites (links provided below).
To learn more about RNASeq in general, please see the links and resources provided on the [**README page**](https://github.com/eisko/RNAseq/blob/master/README.md).


## To get started...
To use the RNASeq tools on your personal laptop, I would recommend using a [**conda environment**](https://docs.conda.io/en/latest/). 

Once conda is installed, you can use the navigator to create a new environment and install the following packages: fastqc, hisat2, stringtie, samtools
Or you can run the following code from the command line:

```
conda create -n RnaSeq
conda activate RnaSeq
conda install fastqc hisat2 stringtie samtools trimmomatic
```
To confirm whether these packages are properly installed run the name of the package as a command in the terminal (e.g. `samtools` or `hisat2`)
If the package is properly installed, you should see a 'help' page for the package.

For example when I run `hisat2` in the terminal I get this as output:
```
No index, query, or output file specified!
HISAT2 version 2.2.0 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
Usage: 
  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <ht2-idx>  Index filename prefix (minus trailing .X.ht2).
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
  <sam>      File for SAM output (default: stdout)

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.

.
.
.
.
.
```

NOTE: These packages can be downloaded on computers using Unix/Linux (e.g. macs), but cannot be downloaded on Windows operating computers.
If you would like to do RNASeq analysis on a windows computer, I would suggest instead using creating a Biowulf account and using Biowulf to do your analysis.

## Step 1: Quality Control
### FastQC
[**FastQC Documentation**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

**Input:** fastq files with sample reads
**Input:** zip and html files containing quality report information

In the terminal, navigate to the folder containing the fastq files you would like to analyze. Run the following command:

`fastqc samplename.fastq`

As output, you should see 2 new files created in the same folder. A `samplename_fastqc.zip` and a `samplename_fastqc.html` file.
Open the html file to view the results in your default web browswer. Using these results, you can inspect the quality of your reads as well as the presence of any adapter sequences.
These results will inform the trimming you may have to perform to improve the quality of your reads.

## Step 2: Trimming
### Trimmomatic
[**Trimmomatic Documentation**](http://www.usadellab.org/cms/?page=trimmomatic)
[**Trimmomatic Manual**](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

**Input:** fastq files
**Output:** fastq files containing trimmed and filtered reads

There is no set way to trim the reads. Each batch of samples may require different trimmed settings.
One example of a trimming strategy:
```
trimmomatic PE -threads 2 \
samplename_R1.fastq samplename_R2.fastq \
samplename_R1_pe.fq samplename_R1_se.fq \
samplename_R1_pe.fq samplename_R1_se.fq \
LEADING:3 \
TRAILING: 3 \
SLIDINGWINDOW:4:20 \
MINLEN:70 \
```
- `PE` indicates the input files are paired end reads
- `-threads 2` indicate to use 2 threads, a setting I use on my mac. In biowulf, this number will depend on the resources allocated, but `$SLURM_CPUS_PER_TASK` can be used.
- `samplename_R1.fastq` and `samplename_R2.fastq` are the forward and reverse reads for a single sample
- `samplename_R1_pe.fq` and `samplename_R1_se.fq` are the designated output files for paired end (`*pe.fq`) and single end (`*se.fq`). For subsequent analysis, I would suggest only using paired end reads.
- `samplename_R2_pe.fq` and `samplename_R2_se.fq` same as above but sequences from reverse instead of forward reads
- `LEADING:3` indicates removing low quality bases from the beginning of the read (bases are N or of quality less than 3)
- `TRAILING: 3` indicates removing low quality bases from the beginning of the read (bases are N or of quality less than 3)
- `SLIDINGWINDOW:4:20` if quality score is less than 20 for any 4 consecutive bases at the beginning or end of read, remove them
- `MINLEN:70` only keep reads that are at least 70 bases long

## Step 3: Alignment
### Hisat2
[**Hisat2 Documentation**](http://daehwankimlab.github.io/hisat2/)

**Input:** fastq file and an indexed reference genome (`.ht2` files) as input
**Output:** Sam files (be warned, these are VERY large files!) and textfile with summary stats

Before running the command, you must have the correct reference genome downloaded. Hisat2 requires an additional step of indexing the reference genome, thus it takes `.ht2` files as input, *not* raw .fa (i.e. fasta files) as input.
You can download the relevant indexed references can be downloaded from the [**Hisat2 webpage**](http://daehwankimlab.github.io/hisat2/download/).
Take note of which reference version you download (e.g. human vs. mouse, GRCm38 vs. mm10) as this can influence your results and requires different annotation files in subsequent steps (i.e. `.gff`/`.gtf` input into **stringtie**)

The following can be run to analyze samples with hisat2:
```
hisat2 --rna-strandness RF --dta -p 2  -x path/to/reference/GRCm38/genome \
-1 path/to/samplename_R1.fq.gz -2 path/to/samplename_R2.fq.gz \
-S path/to/output/folder/samplename.sam \
--summary-file path/to/output/folder/samplename_alignStats.txt
```
- `--rna-strandness RF`
- `--dta`
- `-p`
- `-x path/to/reference/GRCm38/genome`
- `-1 path/to/samplename_R1.fq.gz -2 path/to/samplename_R2.fq.gz`
- `-S path/to/output/folder/samplename.sam`
- `--summary-file path/to/output/folder/samplename_alignStats.txt`


