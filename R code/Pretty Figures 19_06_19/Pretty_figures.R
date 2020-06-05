### Pretty Figures Script
# will output tiff files of pca/volcano - used for poster

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(gplots)
library(ggplot2)
# library(ggbiplot) <- interferes with ballgown functionality?
library(ggfortify)
library(calibrate)
library(cluster)
library(EnhancedVolcano)

setwd("~/code/ballgown2/Pretty Figures")

load("~/code/ballgown2/Vanc_bg.rda")

pData <- read.csv("~/code/ballgown2/Vanc_pdata2.csv")

#filename to use for outputs
file <- "vanc"
# define covariate to make comparisons
# if bg subset = treatment (sham/vanc), then covariate = condition (aff/unaff)
# if bg subset = condition (aff/unaff), then covariate = treatment (sham/vanc)
covariate <- "condition"
# set qvalue used to cut off significant findings
q <- 0.05
fch <- 1
i <- -1 # if sham/vanc group
# i <- 1 # if aff/unaff group
# list of genes of interest (goi) to plot on volcano plot
goi <- c("Mmaa", "Pigr", "Fgf21", "Nnmt")

# colors
# colours <- c("#4370C3", "#70AD46") # for sham c(blue, green)
colours <- c("#4370C3", "#EC7D2F") # for aff c(blue, orange)

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
BG <- subset(bg, "treatment=='vancomycin'", genomesubset=FALSE)
pData <- filter(pData, Treatment=="Vancomycin")
###################################



transcript_fpkm = texpr(BG, 'FPKM')
transcript_fpkm_log2 = log2(transcript_fpkm+1)

pca = prcomp(t(transcript_fpkm_log2))

### see blog for example code
# http://huboqiang.cn/2016/03/03/RscatterPlotPCA
df_out <- as.data.frame(pca$x)
### change labels to male or female
########################### ADJUST for Condition or Treatment
df_out$Treatment <- pData$Treatment
df_out$Sex <- pData$Sex
df_out

# calculate percentages
percentage <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2) # eigen values = square of stddev, traditionally % = proporation of variance, which is based on EIGEN values not std dev

percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )


############ ADJUST GGPLOT AES(COLOR=Condition or Treatment)
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Treatment,shape=Sex))+
  geom_point()+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(colour="grey90", size=0.5),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title = element_text(colour="black", face="bold"),
        axis.text=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"),
        plot.title = element_text(hjust = 0.5))+
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  geom_point(size=5)+
  scale_color_manual(values=colours)

# output tiff for PCA
tiff(filename = paste("PCA_", file, ".tiff"), width = 6, height = 4, units = "in", res = 300)
p
dev.off()



################### VOLCANO PLOT ##########################
################# ADJUST ##########
gene_results <- read.csv(paste("~/code/ballgown2/vanc_gene_results.csv"))

res <- data.frame("ensembl_id" = gene_results$id,
                     "gene_name" = gene_results$gene_name,
                     "log2fc" = i*log2(gene_results$fc), # needed to multiply by -1 to get right up/down orientation in volcano plot
                     "pvalue" = gene_results$pval,
                     "qvalue" = gene_results$qval,
                     "significance" = ifelse(gene_results$qval < q, paste("FDR < ",q), "Not Significant"))

genes <- subset(res, gene_name %in% goi)


######## ADJUST COLOR SCALE?
volcano <- ggplot(res, aes(x = log2fc, y = -log10(pvalue), color=significance)) +
  geom_point(alpha = 0.75) + 
  # scale_color_manual(values = c("red", "black")) + # for sham
  scale_color_manual(values = c("black", "red")) + # for aff
  theme_bw(base_size = 12) + theme(legend.position = "bottom")  +
  geom_point(data=genes, colour="purple") +  # this adds a red point
  geom_label(data=genes, label=genes$gene_name, vjust=1.25, colour="purple")+ # goi text label
  xlab(expression('Log'['2']*' Fold Change'))+
  ylab(expression('-Log'['10']*italic('P')))+
  theme(legend.title=element_blank())+
  xlim(-8,8)+
  ylim(0,6)


# output tiff for volcano
tiff(filename = paste("volcano_", file, ".tiff"), width = 6, height = 4, units = "in", res = 300)
volcano
dev.off()

