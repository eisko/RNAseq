# RNAseq
Scripts and info for RNAseq pipeline

## Scripts
See [**swarm outline**](https://github.com/eisko/RNAseq/blob/master/swarm%20outline) for the pipeline outline and a summary of commands used in swarm scripts on biowulf.

See [**Biowulf_RnaPipeline.sh**](Biowulf_RnaPipeline.sh) for example sbatch script for processing a single sample and recording runtimes.

See [**swarm scripts**](https://github.com/eisko/RNAseq/tree/master/swarm%20scripts) folder for swarm scripts used to run samples on biowulf.
* Swarm scripts were run on sampels named EI01 through EI24

See [**R Code**](https://github.com/eisko/RNAseq/tree/master/R%20code) for R scripts used to process data and make figures.
* [**ballgown_script.R**](R code/ballgown_script.R) - script used to generate differentially expressed genes (DEGs) for comparisons of different conditions
* [**app.R**](R code/app.R) - Interactive shiny app to highlight genes of interest
* [**pretty figures folder**](https://github.com/eisko/RNAseq/tree/master/R%20code/Pretty%20Figures) - files used to generate figures for presentations and posters

---
## Accessory Scripts

See [**folders.sh**](https://github.com/eisko/RNAseq/blob/master/folders.sh) for bash script to generate proper folder structure before running swarm scripts

See [**swarm_files.py**](https://github.com/eisko/RNAseq/blob/master/swarm_files.py) for python script to generate swarm files.
* Script requires you to specify a samples.csv file that conatines 3 columns (Used command line for loop to generate fastq paths):  
   * **SampleIDs** - to label subsequent samples. 
   * **Fastq_Path1** - full path to read1 fastq file. 
   * **Fastq_Path2** - full path to read2 fastq file. 
* also requires several global variables to be defined:  
   * **INDEX** - full path to folder containing hisat indeces - need to run hisat2 build to generate index if one doesn't already exist. 
   * **ANNOTATION** - full path to .gtf file for stringtie to use to assemble/label genes/transcripts. 
   * **PROJECT** - full path to project folder that contains data and proper file structure (see [**swarm outline**](https://github.com/eisko/RNAseq/blob/master/swarm%20outline) for file structure). 

---

## External Links and Resources

Additional resources for learning RNAseq data analysis and necessary tools:
* [**Powerpoint**](https://github.com/eisko/RNAseq/blob/master/RNAseq%20info.pptx) I prepared for lab meeting. Includes overview and some explanation of pipeline
* Overview of RNAseq pipeline: [**Persea et al., 2016**](https://www.nature.com/articles/nprot.2016.095)
* In-depth RNAseq pipeline tutorial: [**Griffith Lab**](https://github.com/griffithlab/rnaseq_tutorial)
* General RNAseq pipeline overview: [**Dave Tang's Blog**](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)
* To implement pipeline on local laptop: [**Conda/bioconda essentials**](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
* Cloud computing on biowulf: [**Biowulf Online Class**](https://hpc.nih.gov/training/intro_biowulf/)
* Interacting with terminal: [**Unix/shell tutorial**](https://www.datacamp.com/courses/introduction-to-shell-for-data-science)

Reviews of RNAseq:
* General overview: [**Wang et al., 2009**](https://www.ncbi.nlm.nih.gov/pubmed/19015660)
* In-depth review: [**Corney, 2013; updated 2018**](https://www.labome.com/method/RNA-seq-Using-Next-Generation-Sequencing.html)

