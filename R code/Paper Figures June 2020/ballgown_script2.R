# ECI 5/18/20
# output gene expression data tables for each condition

# set working directory
setwd("~/code/RnaSeq/ballgown")

# Make sure to necessary packages are installed
library(genefilter)
library(dplyr)
library(devtools)
library(gplots)
library(ggplot2)
library(ggfortify)
library(calibrate)
library(ballgown)

load("~/code/RnaSeq/ballgown/Vanc_bg.rda")

# set working directory
setwd("/Users/iskoec/Documents/quarantine working files/paperfigs pics/rnaseq fig")

#filename to use for outputs
file <- "unaff"

# define covariate to make comparisons
# if bg subset = treatment (sham/vanc), then covariate = condition (aff/unaff)
# if bg subset = condition (aff/unaff), then covariate = treatment (sham/vanc)
covariate <- "treatment"

# set qvalue used to cut off significant findings
q <- 0.1
q2 <- 0.1
fch <- 1

# i = -1 # if sham/vanc
i = 1 # if aff/unaff

# define bg variable for script
# set subset if needed
# treatment==sham ot vancomycin; condition == affected or unaffected (fully spell out)
BG <- subset(bg, "condition=='unaffected'", genomesubset=FALSE)


# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(BG,"rowVars(texpr(BG)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_genes = stattest(bg_filt, feature="gene", covariate=covariate, getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))
results_genes$inverse_fc = (1/results_genes$fc)

# set directory where want to save table
setwd("comparison_results")
# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_genes, paste(file,"_gene_results.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)

