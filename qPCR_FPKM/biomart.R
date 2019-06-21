# practice biomaRt use?

library(biomaRt)


listMarts()
# ensembl=useMart("ensembl")
# datasets <- listDatasets(ensembl)
# head(datasets)
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# or

# want to use: mmusculus_gene_ensembl
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")

# to view possible filters to restrict output for dataset:
filters = listFilters(ensembl)

# possible attributes to retrieve:
attributes = listAttributes(ensembl)


# get ensembl and gene IDs for needed ensembl Ids
ensemblIDs <- read.csv("ensemblIDs.csv", header = FALSE)
ensemblIDs <- as.vector(ensemblIDs$V1)

conversionList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = ensemblIDs,
                        mart = ensembl)

# save conversion list as file, can insert into other R code with load("conversionList.Rda")
save(conversionList, file="conversionList.Rda")


# full conversion list (with every gene)
FullConversionList <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                            mart = ensembl)

#save as object
save(FullConversionList, file = "FullConversionList.Rda")
