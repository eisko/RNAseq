# qPCR process data


import pandas as pd
import numpy as np


data = pd.read_csv("Nnmtall.csv")
hk_gene = "Gapdh" # house keeping gene
gois = ["Nmnt", "Nampt", "Nnmat1"]
pdata = pd.read_csv("Vanc_pdata.csv")

# drop Nan/0 values
data = data.dropna()
data
# average sample values




# group by sample and Primer
sorted <- clean %>% group_by(Sample, Primer)
# calculate means
means <- summarize(sorted, mean(Cq))
# rename column so function doesn't get confused
means <- means %>% rename("avg_Cq" = "mean(Cq)")
# rearrange data so rows = samples and cols = genes
means <- spread(means, Primer, avg_Cq)
