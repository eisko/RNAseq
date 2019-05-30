# RNAseq
Scripts for RNAseq pipeline

See [**swarm outline**](https://github.com/eisko/RNAseq/blob/master/swarm%20outline) for the pipeline outline and a summary of commands used in swarm scripts on biowulf.

See [**Biowulf_RnaPipeline.sh**](Biowulf_RnaPipeline.sh) for example sbatch script for processing a single sample and recording runtimes.

See [**swarm scripts**](https://github.com/eisko/RNAseq/tree/master/swarm%20scripts) folder for swarm scripts used to run samples on biowulf.
* Swarm scripts were run on sampels named EI01 through EI24

See [**swarm_files.py**](https://github.com/eisko/RNAseq/blob/master/swarm_files.py) for python script to generate swarm files.
* Script requires you to specify a samples.csv file that conatines 3 columns (Used command line for loop to generate fastq paths):  
   * **SampleIDs** - to label subsequent samples. 
   * **Fastq_Path1** - full path to read1 fastq file. 
   * **Fastq_Path2** - full path to read2 fastq file. 
* also requires several global variables to be defined:  
   * **INDEX** - full path to folder containing hisat indeces - need to run hisat2 build to generate index if one doesn't already exist. 
   * **ANNOTATION** - full path to .gtf file for stringtie to use to assemble/label genes/transcripts. 
   * **PROJECT** - full path to project folder that contains data and proper file structure (see [**swarm outline**](https://github.com/eisko/RNAseq/blob/master/swarm%20outline) for file structure). 

See [**R Code**](https://github.com/eisko/RNAseq/tree/master/R%20code) for R scripts used to process data (ballgown_script.R) and make figures (need to include more).
* **ballgown_script.R** - script used to generate DEGs for comparisons of different conditions
* **app.R** - Interactive shiny app to highlight genes of interest
* pretty figures ... - files used to generate figures for presentations and posters

