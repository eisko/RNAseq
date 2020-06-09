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

NOTE: I would recommend at this step using a simple `samplename` prefix and numbering sample 01, 02, 03,...,10,11,12... to make subsequent ballgown data analysis  easier.

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
[**Ballgown Vignette**](https://www.bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html) <- highly recommend checking this out. Great vignette and info on how to use the package and what else it can do! \
[**Ballgown Paper**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4792117/) 

Below, you can find the code I used to analyze the data. I have annotated it, but it may need tweaking to fit your data and your analysis.

Load the following packages:
```r
library("genefilter")
library("dplyr")
library("devtools")
library("ballgown")
```

Use the output files of stringtie to create the ballgown object. The output for stringtie should look something like this:
```
experiment_folder/
    sample01/
        e2t.ctab
        e_data.ctab
        i2t.ctab
        i_data.ctab
        t_data.ctab
    sample02/
        e2t.ctab
        e_data.ctab
        i2t.ctab
        i_data.ctab
        t_data.ctab
    ...
    sample20/
        e2t.ctab
        e_data.ctab
        i2t.ctab
        i_data.ctab
        t_data.ctab
```
Where experiment_folder contains a labeled folder for each sample processed, and each sample folder (e.g. `sample01`) contain 5 `.ctab` files. Additional files created by stringtie in these folders (e.g. a `.gtf` or `.tsv` files) will not interfere with the generation of the ballgown object.

If the data is organized in that way, generating the ballgown object becomes simple:
```r
bg <- ballgown(dataDir="path/to/experiment_folder", samplePattern='sample', meas='all')
```
- `dataDir=` provide a string of the path leading and including the `experiment_folder`
- `samplePattern=` provide the prefix used for naming samples (e.g. I used `EI_01` to `EI_24`, so I supplied `samplePattern='EI_'`)
- `meas=` this indicates the expression measure used for creating the ballgown object. Other objects include `'FPKM'`, `'cov'`, and more.

If your samples are not numerically or simply labelled, you can instead create a vector containing the paths to all your samples and use this to generate your ballgown object:
```r
sample_list <- c("/path/to/experiment_folder/sample01",
                 "/path/to/experiment_folder/sample02",
                 "/path/to/experiment_folder/sample03",
                 "/path/to/experiment_folder/sample04",
                 .
                 .
                 .
                 .
                 )
bg <- ballgown(samples=sample_list, meas='all')
```


This step will likely take awhile and create a large object. I suggest saving it as a RDS object, so you don't have to repeated build the ballgown object
```r
saveRDS("path/to/ballgown/output/bg.rds")
```

Now that the ballgown object is created, we can inspect the data by generating a table of transcript or gene expression values. See the [**Ballgown Vignette**](https://www.bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html) for more information and options.\
- `texpr(bg)` will generate a table where columns are individual samples, rows are transcripts, and each entry is the FPKM observed in each sample for that specified transcript
- `gexpr(bg)` will generate a similar table but each row is a gene instead of a transcript
Example:
```r
texpr(bg)
```

```r
 FPKM.EI_01 FPKM.EI_02    . . . FPKM.EI_10 FPKM.EI_11 FPKM.EI_12 FPKM.EI_13 FPKM.EI_14
1          0          0   . . .   0.000000          0          0          0   0.000000
2          0          0   . . .   0.000000          0          0          0   0.000000
3          0          0   . . .   0.000000          0          0          0   0.018331
4          0          0   . . .   0.000000          0          0          0   0.000000 
5          0          0   . . .   0.014456          0          0          0   0.000000
6          0          0   . . .   0.000000          0          0          0   0.000000
```

Notice how there are a lot of 0s, this will be addressed in the filtering step below.


The next step is to import metadata (i.e. 'phenotypic data') and assign it to the ballgown object. This information is generally stored in a `.csv` or a `.tsv` file, which can be created in excel. My metadata file looks something like this:
```
sample_id,condition,treatment,sex,group
EI_01,affected,sham,F,Af_Sham
EI_02,affected,sham,F,Af_Sham
EI_03,affected,sham,F,Af_Sham
EI_04,affected,sham,M,Af_Sham
EI_05,affected,sham,M,Af_Sham
EI_06,affected,sham,M,Af_Sham
EI_07,affected,vancomycin,F,Af_Vanc
EI_08,affected,vancomycin,F,Af_Vanc
EI_09,affected,vancomycin,F,Af_Vanc
EI_10,affected,vancomycin,M,Af_Vanc
EI_11,affected,vancomycin,M,Af_Vanc
EI_12,affected,vancomycin,M,Af_Vanc
EI_13,unaffected,sham,F,Un_Sham
EI_14,unaffected,sham,F,Un_Sham
EI_15,unaffected,sham,F,Un_Sham
EI_16,unaffected,sham,M,Un_Sham
EI_17,unaffected,sham,M,Un_Sham
EI_18,unaffected,sham,M,Un_Sham
EI_19,unaffected,vancomycin,F,Un_Vanc
EI_20,unaffected,vancomycin,F,Un_Vanc
EI_21,unaffected,vancomycin,F,Un_Vanc
EI_22,unaffected,vancomycin,M,Un_Vanc
EI_23,unaffected,vancomycin,M,Un_Vanc
EI_24,unaffected,vancomycin,M,Un_Vanc
```

NOTE: make sure the first column (`sample_id`) *exactly* matches the sample names used to label the folders containing ballgown input data (e.g. `sample01` in the example file structure above).

Import the metadata and assign it to the pData slot of the ballgown object.
```r
metadata <- read.csv("path/to/metadata.csv")

pData(bg) <- metadata
pData(bg)
```

Determine what samples you want to compare and across what variable. For example, if I wanted to determine what genes are differentially expressed in the disease state, I would need to restrict my analysis to samples under the `sham` treatment and then compare the `affected` and `unaffected` conditions.

To restrict the analysis to only animals under the sham treatment, I would create a subsetted ballgown object:
```r
bg_subset <- subset(bg, "treatment=='sham'", genomesubset=FALSE)
```

Before completing the statistical analysis, low abundance transcripts must be filtered out to reduce noise and the sparseness in our dataset. The following code will discard genes in which the variance across samples is less than 1 (i.e. the gene has very low abundance in the overally dataset):
```r
bg_filt = subset(bg_subset,"rowVars(texpr(bg_subset)) > 1", genomesubset=TRUE)
```

If you reinspect the ballgown object you will notice many less 0s (i.e. reducing the 'sparseness' of the dataset):
```r
texpr(bg_filt)
```
```
FPKM.EI_01 FPKM.EI_02 FPKM.EI_03 FPKM.EI_04 FPKM.EI_05 FPKM.EI_06 FPKM.EI_13 FPKM.EI_14 FPKM.EI_15
46   5.430120   4.722524   1.880314   3.575279   4.019314   7.252697   2.528504   2.482140   3.035233
47   7.452411   8.411320   5.454850   7.409715   9.295838   5.133492   5.535621   5.947452   6.402004
48   2.800960   5.133840   1.824513   4.669214   6.825549   6.026739   6.467938   3.099066   5.983899
52   1.927685   8.972164   3.501288   3.722792   5.407292   0.000000   4.631611   4.008284   6.711938
55  38.717491  55.873985  26.445629  47.138226  43.783592  56.063797  47.812008  60.119213  56.030407
63  26.751963  30.930845  53.045113  25.794495  22.391626  24.219151  17.764759  16.974792  19.402292
```



Next, we will compare expression data across a variable in the metadata. In my case, I have restricted to the samples to `sham` treatment, and I would like to compare across the `'condition'` (i.e. `affected` vs. `unaffected`). Ballgown uses the normalization and anlysis implemented in the `limma` package. To learn more and find more references on the statistical analysis used in ballgown see the [**Ballgown Paper**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4792117/) or [**Ballgown Vignette**](https://www.bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html). To generate the statistical analysis run:
```r
results_genes = stattest(bg_filt, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
```
- `bg_filt` the ballgown object to use for the analysis
- `feature="gene"` the feature to be used for the analysis. Other options include `"transcript"` or `"exon"`
- `covariate="condition"` a column specified in the metadata/pData to sort and compare the samples into different groups. There must be at least 2 attributes specified in the variable (in my case `affected` and `unaffected`) in order for ballgown to compare across 2 different groups.
- `getFC=TRUE` ensures fold change values are returned
- `meas="FPKM"` specifies what measure to use for the analysis, in this case `FPKM` is used.

The output of `stattest()` is a table containing the results of the statistical test. It should look something like this:
```r
head(results_genes)
```

```r
 feature                 id        fc        pval      qval
1    gene ENSMUSG00000000001 1.3316766 0.051766131 0.3137894
2    gene ENSMUSG00000000031 0.3162352 0.343807799 0.6383109
3    gene ENSMUSG00000000049 0.7694331 0.223577924 0.5329133
4    gene ENSMUSG00000000056 1.2834758 0.187591088 0.5020498
5    gene ENSMUSG00000000078 1.5242118 0.404437300 0.6851108
6    gene ENSMUSG00000000088 0.5033407 0.007286408 0.1807343
```

To append gene names to the results table:
```r
# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')

# Extract only gene name and ensembl ID columns
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# add gene names to results table based on matching Ensembl ID
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

head(results_genes)
```

```r
                  id feature        fc        pval      qval gene_name
1 ENSMUSG00000000001    gene 1.3316766 0.051766131 0.3137894     Gnai3
2 ENSMUSG00000000031    gene 0.3162352 0.343807799 0.6383109       H19
3 ENSMUSG00000000049    gene 0.7694331 0.223577924 0.5329133      Apoh
4 ENSMUSG00000000056    gene 1.2834758 0.187591088 0.5020498      Narf
5 ENSMUSG00000000078    gene 1.5242118 0.404437300 0.6851108      Klf6
6 ENSMUSG00000000088    gene 0.5033407 0.007286408 0.1807343     Cox5a
```

Make sure to save the results as a `.csv` or `.tsv` file. This file can be opened and inspected in excel and can be used as supplementary data upon publication.
```r
write.table(results_genes, "path/to/output/filename_gene_results.csv", sep=",", quote=FALSE, row.names = FALSE)
```
