# RNAseq
Scripts and info for RNAseq pipeline

## Tutorial
See [**Emily's RnaSeq Tutorial**](https://github.com/eisko/RNAseq/blob/master/Emily's%20RNASeq%20Tutorial.md) for a tutorial I generated myself. It is not completely exhaustive and does not go in depth in explaining what each program does and how. If you would like to know more about what each program does and how, I suggest visiting the software page and looking at the documentation or looking up the published article associated with each software tool.

## Scripts
See [**swarm outline**](https://github.com/eisko/RNAseq/blob/master/swarm%20outline) for the pipeline outline and a summary of commands used in swarm scripts on biowulf.

See [**Biowulf_RnaPipeline.sh**](Biowulf_RnaPipeline.sh) for example sbatch script for processing a single sample and recording runtimes.

See [**swarm scripts**](https://github.com/eisko/RNAseq/tree/master/swarm%20scripts) folder for swarm scripts used to run samples on biowulf.
* Swarm scripts were run on sampels named EI01 through EI24

See [**R Code**](https://github.com/eisko/RNAseq/tree/master/R%20code) for R scripts used to process data and make figures.
* [**Paper figures**](https://github.com/eisko/RNAseq/tree/master/R%20code/Pretty%20Figures) - Contains R Scripts and markdown files used to generate data and figures for paper
* [**app.R**](https://github.com/eisko/RNAseq/blob/master/R%20code/app.R) - Interactive shiny app to highlight genes of interest
* For R tutorial with ballgown see scripts copied from [**Griffith Tutorial**](https://github.com/griffithlab/rnaseq_tutorial/tree/master/scripts) or find relevant scripts in [**Griffith R Scripts**]() folder.


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

## More Links and Resources

Great overview of bioinformatics basics:
* [**iab online textbook**](https://readiab.org/introduction.html) is a great interactive, free resources for learning about bioinformatics and how to manipulate sequences in python

Additional resources for learning RNAseq data analysis and necessary tools:
* [**Powerpoint**](https://github.com/eisko/RNAseq/blob/master/RNAseq%20info.pptx) I prepared for lab meeting. Includes overview and some explanation of pipeline
* Overview of RNAseq pipeline: [**Pertea et al., 2016**](https://www.nature.com/articles/nprot.2016.095)
* In-depth RNAseq pipeline tutorial: [**Griffith Lab**](https://github.com/griffithlab/rnaseq_tutorial)
  * UPDATE: new/updated tutorial from Griffith lab that reviews mroe complex tools that goes beyond the Pertea pipeline: [**RNAseq**](https://rnabio.org/course/)
* General RNAseq pipeline overview: [**Dave Tang's Blog**](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)
* To implement pipeline on local laptop: [**Conda/bioconda essentials**](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html)
* Cloud computing on biowulf: [**Biowulf Online Class**](https://hpc.nih.gov/training/intro_biowulf/)
* Interacting with terminal: [**Unix/shell tutorial**](https://www.datacamp.com/courses/introduction-to-shell-for-data-science)
* Another great teaching resource for learning theory and implementation of RNAseq: [**Harvard Chan Bioinformatics Core Training Program**](https://hbctraining.github.io/main/#introduction-to-next-generation-sequencing-ngs-analysis-series)

Reviews of RNAseq:
* General overview: [**Wang et al., 2009**](https://www.ncbi.nlm.nih.gov/pubmed/19015660)
* In-depth review: [**Corney, 2013; updated 2018**](https://www.researchgate.net/publication/306291842_RNA-seq_Using_Next_Generation_Sequencing_A_comprehensive_review_of_RNA-seq_methodologies)

