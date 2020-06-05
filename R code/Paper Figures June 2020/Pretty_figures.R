### Pretty Figures Script
# will output tiff files of pca/volcano - used for poster
# last update: ECI 5/14/20

library(genefilter)
library(dplyr)
library(devtools)
library(gplots)
library(ggplot2)
library(ggfortify)
library(calibrate)
library(cluster)
library(ggrepel)
library(ballgown)


setwd("~/Documents/quarantine working files/paperfigs pics/rnaseq fig")

load("~/code/RnaSeq/ballgown/Vanc_bg.rda")

pData <- read.csv("~/code/RnaSeq/ballgown/Vanc_pdata2.csv")


#filename to use for outputs
################# ADJUST
file <- "unaff"

# define covariate to make comparisons
########################### ADJUST
# covariate <- "condition" # if sham/vanc group
covariate <- "treatment" # if aff/unaff

# set qvalue used to cut off significant findings
q <- 0.1
fch <- 1
##############################ADJUST
# i <- -1 # if sham/vanc
i <- 1 # if aff/unaff group

# list of genes of interest (goi) to plot on volcano plot
goi <- c("Mmaa", "Fgf21")

# colors
blue <- "#4370C3"
green <- "#70AD46"
orange <- "#EC7D2F"
yellow <- "#FDB40B"

################################# ADJUST
# colours <- c(blue, green) # for sham
# colours <- c(blue, orange) # for aff
# colours <- c(orange, yellow) # for vanc
# colours <- c(green, yellow) # for unaff

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             # panel.grid.major = element_line(colour="grey90", size=0.5),
             # panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title = element_text(colour="black", face="bold"),
             axis.text=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"),
             plot.title = element_text(hjust = 0.5))




# ADJUST FOR SUBSET#######################
#### fully spell out, treatment==sham or vancomycin; condition == affected or unaffected
BG <- subset(bg, "condition=='unaffected'", genomesubset=FALSE)
pData <- filter(pData, Condition=="Unaffected")
###################################

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(BG,"rowVars(texpr(BG)) > 1", genomesubset=TRUE)

bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

results_genes = stattest(bg_filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# write.csv(results_genes, file=paste0("./Volcano tables/", file, "_gene_results.csv"))


gene_results <- results_genes

################### VOLCANO PLOT ##########################
res <- data.frame("ensembl_id" = gene_results$id,
                     "gene_name" = gene_results$gene_name,
                     "log2fc" = i*log2(gene_results$fc), # needed to multiply by -1 to get right up/down orientation in volcano plot
                     "pvalue" = gene_results$pval,
                     "qvalue" = gene_results$qval,
                     "significance" = ifelse(gene_results$qval < q, paste("FDR < ",q), "Not Significant"))

# create lists of genes to highlight in volcano graph
genes <- subset(res, gene_name %in% goi)
res2 <- arrange(res, pvalue)
genes_15 <- res2[1:15,]


######## ADJUST COLOR SCALE
volcano <- ggplot(res, aes(x = log2fc, y = -log10(pvalue), color=significance)) +
  geom_point(alpha = 0.75) + 
  ################ ADJUST:
  scale_color_manual(values = c("red", "black")) + # for sham/vanc/unaff
  # scale_color_manual(values = c("black", "red")) + # for aff
  theme_bw(base_size = 8) + theme(legend.position = "bottom")  +
  geom_label_repel(data=genes_15, label=genes_15$gene_name, colour="black")+ 
  geom_point(data=genes, colour="purple") +  # this adds a red point
  geom_label_repel(data=genes, label=genes$gene_name, colour="purple")+ # goi text label
  xlab(expression('Log'['2']*' Fold Change'))+
  ylab(expression('-Log'['10']*italic('P')))+
  theme(legend.title=element_blank())+
  xlim(-8,8)+
  ylim(0,6) +
  theme(axis.title = element_text(face="bold", size=20),
        legend.text = element_text(size=20))


# output tiff for volcano
tiff(filename = paste("pics/volcano_", file, ".tiff"), width = 6, height = 6, units = "in", res = 300)
volcano
dev.off()

# transcript_fpkm = texpr(bg_filt, 'FPKM')
# transcript_fpkm_log2 = log2(transcript_fpkm+1)
# 
# pca = prcomp(t(transcript_fpkm_log2))
# 
# ### see blog for example code
# # http://huboqiang.cn/2016/03/03/RscatterPlotPCA
# df_out <- as.data.frame(pca$x)
# ### change labels to male or female
# ########################### ADJUST for Condition (sham/vanc) or Treatment (aff/unaff)
# df_out$Treatment <- pData$Treatment
# df_out$Sex <- pData$Sex
# df_out
# 
# # calculate percentages
# percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2) # eigen values = square of stddev, traditionally % = proporation of variance, which is based on EIGEN values not std dev
# 
# percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
# 
# 
# ############ ADJUST GGPLOT AES(COLOR=Condition (sham/vanc) or Treatment (aff/unaff))
# p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Treatment,shape=Sex))+
#   geom_point()+
#   theme(panel.background = element_blank(),
#         panel.border=element_rect(fill=NA),
#         panel.grid.major = element_line(colour="grey90", size=0.5),
#         panel.grid.minor = element_blank(),
#         strip.background=element_blank(),
#         axis.title = element_text(colour="black", face="bold"),
#         axis.text=element_text(colour="black"),
#         axis.ticks=element_line(colour="black"),
#         plot.margin=unit(c(1,1,1,1),"line"),
#         plot.title = element_text(hjust = 0.5))+
#   xlab(percentage[1]) + 
#   ylab(percentage[2])+
#   geom_point(size=5)+
#   scale_color_manual(values=colours)
# 
# # output tiff for PCA
# tiff(filename = paste("pics/PCA_", file, ".tiff"), width = 6, height = 4, units = "in", res = 300)
# p
# dev.off()
