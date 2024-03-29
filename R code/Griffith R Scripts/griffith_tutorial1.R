#Jason Walker, jason.walker[AT]wustl.edu
#Malachi Griffith, mgriffit[AT]wustl.edu
#Obi Griffith, obigriffith[AT]wustl.edu
#The Genome McDonnell Institute, Washington University School of Medicine

#R tutorial for Informatics for RNA-sequence Analysis workshops

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

#load bg object
load("Vanc_bg.rda")

# Load phenotype data from a file we saved in the current working directory
pheno_data = read.csv("Vanc_pdata.csv")

pData(bg) = pheno_data

# Load ballgown data structure and save it to a variable "bg"
#bg <- ballgown(dataDir = "/data/iskoec/Vanc_RnaSeq/stringtie", samplePattern = 'EI_', pData=pheno_data, meas='all')

# Display a description of this object
bg


# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])


# Perform differential expression (DE) analysis on TREATMENT (antibiotic vs. sham) with no filtering
# NOTE: coveriate variable determines comparison (can set to "condition", "treatment", or "sex")
results_transcripts = stattest(bg, feature="transcript", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "Mutant_vs_Het_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "Mutant_vs_Het_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "Mutant_vs_Het_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "Mutant_vs_Het_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts, "Mutant_vs_Het_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "Mutant_vs_Het_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Exit the R session
quit(save="no")