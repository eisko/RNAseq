# process Nnmt values?

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/code/RnaSeq/ballgown/Nnmt")
nnmtData <- read.csv("Nnmt_FPKM.csv")
nnmtDataT <- t(nnmtData)
nnmtDF <- as.data.frame(nnmtDataT[2:25,])
nnmtDF_tbl <- tibble::rownames_to_column(nnmtDF, var = "sample_id")
names(nnmtDF_tbl) <- c("sample_id", "Nnmt_FPKM")


pData <- read.csv("~/code/RnaSeq/ballgown/Vanc_pdata.csv")



nnmtDatafull <- left_join(nnmtDF_tbl, pData, by = c("sample_id"= "sample_id"))

nnmtDatafull <- arrange(nnmtDatafull, treatment)

write.csv(nnmtDatafull, file="nnmt_fpkm_arranged.csv")


# p <- ggplot(nnmtDatafull, aes(x= "sample_id", y= "Nnmt_FPKM"))+
#   geom_bar()
# p
