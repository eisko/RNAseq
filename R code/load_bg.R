library(ballgown)

bg <- ballgown(dataDir = "~/code/RnaSeq/stringtie", samplePattern = 'EI_', meas='all')

pData(bg) = read.csv("Vanc_pdata.csv")

save(bg, file='Vanc_bg.rda')
