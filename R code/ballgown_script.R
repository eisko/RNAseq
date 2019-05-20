# R script for processing DE Rna Seq data, or subset of data
# using ballgown package
# 2/12/19 ECI

# NOTE: Edit code so that it takes input: subset = " ", covariate = " "


# set working directory
setwd("~/code/RnaSeq/ballgown")

# Make sure to necessary packages are installed
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(gplots)
library(ggplot2)
library(ggbiplot)
library(ggfortify)
library(calibrate)
#library(pca3d) # <-- package not available for mac os


#load bg object
# see "load_bg.R" for making of ballgown (bg) object and subsets of bg object
load("Vanc_bg.rda")

# set working directory
setwd("~/code/RnaSeq/ballgown")

# define bg variable for script
# set subset if needed
BG <- subset(bg, "treatment=='sham'", genomesubset=FALSE)

#filename to use for outputs
file <- "sham"

# define covariate to make comparisons
# if bg subset = treatment (sham/vanc), then covariate = condition (aff/unaff)
# if bg subset = condition (aff/unaff), then covariate = treatment (sham/vanc)
covariate <- "condition"

# set qvalue used to cut off significant findings
q <- 0.1
q2 <- 0.1
fch <- 1

# if sham/vanc group, set to -1
# if aff/unaff group, set to 1
i = -1


# list of genes of interest (goi) to plot on volcano plot
goi <- c("Mmaa", "Pigr", "Fgf21", "NMNT")


# set up output folder
dir.create(file, recursive = TRUE)
setwd(file)


#extract transcript FPKM values
transcript_fpkm = texpr(BG, 'FPKM')
transcript_fpkm_log2 = log2(transcript_fpkm+1) #+1 prevents nonsense values for no gene expression


# BOXPLOT
jpeg(filename= paste(file, "_boxplot.jpeg", sep=""), width = 676, height = 473)
boxplot (transcript_fpkm_log2, col=as.numeric(pData(BG)$group),las=2,ylab='log2(FPKM+1)')
dev.off()

#PCA
jpeg(filename= paste(file, "_pca.jpeg", sep=""), width = 676, height = 473)
pca = prcomp(t(transcript_fpkm_log2))
autoplot(pca,data=pData(BG),colour=covariate, label=T, label.label="sample_id") + ggtitle(paste(file,"PCA"))
# ggbiplot(pca, ellipse = T )
dev.off()

## write out tables of gene expression

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(BG,"rowVars(texpr(BG)) > 1", genomesubset=TRUE)

# Generate .gtc file for GSEA application (gene expression output)
gene_expr <- gexpr(bg_filt)
gene_expr <- cbind(NA, gene_expr)
header <- data.frame(c("#1.2", nrow(gene_expr)), c(NA, ncol(gene_expr)))

write.table(header, paste(file,"_gexpr.tsv",sep=""), sep="\t", quote=FALSE)
write.table(gene_expr, paste(file,"_gexpr.tsv",sep=""), sep="\t", quote=FALSE, append = TRUE)
write.table(pData(bg_filt), paste(file,"_pData.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)





# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_genes = stattest(bg_filt, feature="gene", covariate=covariate, getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))
results_genes$inverse_fc = (1/results_genes$fc)

# capture all up and down regulated genes
genes_up = i*log2(results_genes[,"fc"]) > 0
genes_up = results_genes[genes_up,]
genes_down = i*log2(results_genes[,"fc"]) < 0
genes_down = results_genes[genes_down,]

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_genes, paste(file,"_gene_results.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)

# Output table of either up or down regulated genes
write.table(genes_up, paste(file,"_up_gene_results.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)
write.table(genes_down, paste(file,"_down_gene_results.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)


# VOLCANO
# Make a basic volcano plot
# NOTE: need to can input -log2(fc) or log2(fc), depending on stattest results
# stattest function will calculate fc with second_level/first_level
jpeg(filename= paste(file, "_volcano.jpeg", sep=""), width = 676, height = 473)

with(results_genes, plot(i*log2(fc), -log10(pval), pch=20, main=paste(file, "Volcano plot"), xlim=c(-10,10), ylim=c(0,11))) 

# Add colored points: orange if q value<q, green if log2FC>1.5, red if log2FC<-1.5, blue if both)
with(subset(results_genes, qval<q ), points(i*log2(fc), -log10(pval), pch=20, col="orange"))
with(subset(results_genes, i*log2(fc)> fch), points(i*log2(fc), -log10(pval), pch=20, col="green"))
with(subset(results_genes, i*log2(fc)< -fch), points(i*log2(fc), -log10(pval), pch=20, col="red"))
with(subset(results_genes, qval<q & abs(log2(fc))>fch), points(i*log2(fc), -log10(pval), pch=20, col="blue"))

# plot genes of interest (goi) with one command
goi <- c("Nnmt", "Mmaa", "Fgf21", "Pigr", "Tuba8")
goi_val <- filter(results_genes_sham, gene_name %in% goi)
with(goi_val, points(-log2(fc), -log10(pval), pch=20, col="purple"))
with(goi_val, textxy(-log2(fc), -log10(pval), labs=gene_name,col="purple",cex=.8))

# Add text saying how many up, down, and total DE genes were found
text(7, y=9, labels = paste("Up DE genes: ", nrow(genes_up)), cex = 0.9)
text(-7, y=9, labels = paste("Down DE genes: ", nrow(genes_down)), cex = 0.9)
text(0, y=8, labels = paste("Total DE genes: ", nrow(results_genes)), cex = 0.9)

# Add labels
with(subset(results_genes, qval<.05 & abs(log2(fc))>1), textxy(-log2(fc), -log10(pval), labs=gene_name, cex=.8))

dev.off()


# # write out tables of DE genes w/ fold change cutoff of fch
# subset_up = subset(results_genes, qval<q2 & i*log2(fc) > fch)
# write.table(subset_up, paste(file,"_fc",fch,"_q",q,"_up_gene_results_filtered.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)
# 
# subset_down = subset(results_genes, qval<q2 & i*log2(fc)< -fch)
# write.table(subset_down, paste(file,"_fc",fch,"_q",q,"_down_gene_results_filtered.tsv",sep=""), sep="\t", quote=FALSE, row.names = FALSE)



