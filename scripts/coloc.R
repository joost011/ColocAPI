library(jsonlite)
library(coloc)
library(rjson)

getwd()

args <- commandArgs(trailingOnly=TRUE)
gene <- args[1]
file_id <- args[2]

data <- fromJSON(file = paste("storage/processed_files/", file_id, "/", gene, ".json", sep=''))

data$gwas$type <- 'cc'
data$gwas$N <- 276020
data$gwas$s <- 0.1284

data$eqtls$type <- 'quant'

# Perform colocalization
res <- coloc.abf(dataset1=data$gwas, dataset2=data$eqtls)

datasets <- list(datasets = list(location = data$location, gwas = data$gwas, eqtls = data$eqtls, coloc = res$results, posterior = res$summary))

if (res$summary[6] >= 0.8){
  write(toJSON(datasets), paste("storage/out_files/", file_id, "/", gene, ".json", sep=''))
}