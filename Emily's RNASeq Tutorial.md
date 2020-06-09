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

**Input:** fastq files with sample reads \
**Output:** zip and html files containing quality report information

In the terminal, navigate to the folder containing the fastq files you would like to analyze. Run the following command:

`fastqc samplename.fastq`

As output, you should see 2 new files created in the same folder. A `samplename_fastqc.zip` and a `samplename_fastqc.html` file.
Open the html file to view the results in your default web browswer. Using these results, you can inspect the quality of your reads as well as the presence of any adapter sequences.
These results will inform the trimming you may have to perform to improve the quality of your reads.

## Step 2: Trimming
### Trimmomatic
[**Trimmomatic Documentation**](http://www.usadellab.org/cms/?page=trimmomatic) \
[**Trimmomatic Manual**](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

**Input:** fastq files \
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

**Input:** fastq file and an indexed reference genome (`.ht2` files) as input \
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
- `--rna-strandness RF` indicates reads are paired whether the first file (`-1`) is forward or reverse read
- `--dta` report alignments tailored for downstream assemblers. Needed to ensure output can be used for **stringtie**
- `-p` indicates number of threads/cpus to use to process data. Home laptop I use 2, biowulf I use $SLURM_CPUS_PER_TASK
- `-x path/to/reference/GRCm38/genome` details path to reference genome. Provide full path to where genome.ht2 files are located. Make sure to also specify the shared prefix of all the reference files in the folder, in this case the prefix of each file is `genome`. I recommend providing the FULL path
- `-1 path/to/samplename_R1.fq.gz -2 path/to/samplename_R2.fq.gz` provide full path to forward and reverse reads for a single sample
- `-S path/to/output/folder/samplename.sam` provide the path and name of file you want as output. If not specified, hisat will create a `samplename.sam` file in the directory where you run the command
- `--summary-file path/to/output/folder/samplename_alignStats.txt` will create a textfile with a summary of the read alignments (i.e. how many reads did or did not align)

## Step 4: Sort
### Samtools
[**Samtools Documentation**](http://www.htslib.org/) \
[**Samtools Manual**](http://www.htslib.org/doc/samtools.html)

**Input:** sam file as input (be warned, these are VERY large files!) \
**Output:** bam file (a compressed version of the sam file)

The purpose of this step is to compress the Sequence Alignment/Map (SAM) file into a bam (binary version of a sam file) format. Bam files are much smaller and easier for the computer to handle. Each line in the file documents a read alignment to the reference genome. You can learn morea bout sam/bam file format [**here**](https://samtools.github.io/hts-specs/SAMv1.pdf).

This step also sorts the alignemnts by where the order they appear in the genome (alignments to chromosome 1 first, chr2...). This step is needed for better `stringtie` assembly.
```
samtools sort -@ 2 -o /path/to/output/samplename_sorted.bam /path/to/input/samplename.sam 
```
- `sort` the function you want samtools to do. Samtools can perform other useful function (like `samtools view`) refer to the manual to explore other samtools functions that may be useful
- `-@` specifies how many threads to use to run code
- `-o /path/to/output/samplename_sorted.bam` specifies the output `.bam` file to be created
- `/path/to/input/samplename.sam` specify the input `.sam` file to perform the function on. Note: the input does not require a flag

After creating the `.bam` files, I would recommend deleting (using the unix `rm` function) the `.sam` files as they take up a lot of space and all subsequent steps can be run using the `.bam` files. 

If you would like to visually inspect the `.bam` files using [**IGV**](http://software.broadinstitute.org/software/igv/), you first need to create a bam index using the samtools index function (`samtools samtools index /path/to/samplename_sorted.bam`).

## Step 5: Assembly
### Stringtie
[**Stringtie Documentation**](https://ccb.jhu.edu/software/stringtie/) \
[**Stringtie Manual**](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

**Input:** bam file and genome annotation file (`.gff` or `.gtf` file will work) \
**Output:** several files in a `samplename` folder
- `sample_transcripts.gtf` contains information on all of the transcripts identified using stringtie
- `.ctab` files (there should be 5 of them per file). These files are used as input for ballgown
- `samplename_gene_abundances.tsv` a matrix where columns are samples, rows are genes, and entries are the raw counts of reads/gene/sample. This can be used as an input count matrix into RNASeq data processing packages in R

The purpose of this step is to 'count' the reads per intron/exon/transcript/gene. 

```
stringtie -p 2 -G \path\to\genome_annotation.gtf -e -B \
-o samplename_transcripts.gtf \
-A samplename_gene_abundances.tsv \
path/to/input/samplename_sorted.bam
```
- `-p 2` number of threads
- `-G \path\to\genome_annotation.gtf` specifies that a genome assembly file should be used to do the assembly. If no annotation files is given and `-G` is not specified, stringtie will do a 'de novo' assembly, i.e. will assemble transcripts without a reference.
- `-e` Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option. This is recommended by the makers of Stringtie for downstream ballgown/edgeR/DEseq analysis.
- `B` specifies to create files (i.e. `.ctab` files) needed for ballgown analysis
- `-o samplename_transcripts.gtf` specifies output `.gtf` file containing info on transcripts identified in the data, does not contain information on 'count' data i.e. how many reads found per sample

**************************
# R Code Tutorial

In R, you will be able to 1) create a results table of differential expression, 2) create a volcano plot, 3) create a PCA plot.

The main packages used to analyze the data include `ballgown` and `edgeR`. Other packages used to analyze RnaSeq data include `DESeq2`, but I have not explored and will not be using the package. Each of these packages use slightly different statistical models and methods to analyze the data, and it may be worth doing some brief reserach to determine which package/analysis method best suits your data.

## Before you start...

Many of these packages (e.g. `ballgown` and `edgeR`) can be installed through [**bioconductor**](https://www.bioconductor.org/). I would suggest *not* installing these packages through conda as this can cause glitches to occur. Instead, download Rstudio to your laptop and use bioconductor (`BiocManager::install("package")`) in the R console to download packages. Using bioconductor rather than the base `install()` function will ensure the packages and all its dependencies are configured corrrectly.

## Step 1: Differential Expression analysis
### Ballgown
[**Ballgown Installation**](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html) \
[**Ballgown Vignette**](https://www.bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html) <- highly recommend checking this out. Great vignette and info on how to use the package and what else it can do!

Below, you can find the code I used to analyze the data. I have annotated it, but it may need tweaking to fit your data and your analysis.

Load the following packages:
```r
library(genefilter)
library(dplyr)
library(devtools)
library(ballgown)
```




