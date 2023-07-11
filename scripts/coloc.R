library(jsonlite)
library(coloc)
library(rjson)

getwd()

args <- commandArgs(trailingOnly=TRUE)
gene <- args[1]
file_id <- args[2]
coloc_type <- args[3]

data$gwas$type <- coloc_type
data$gwas$N <- as.integer(args[4])

if (coloc_type == 'cc'){
  data$gwas$s <- as.numeric(args[5])
} 

data <- fromJSON(file = paste("storage/processed_files/", file_id, "/", gene, ".json", sep=''))

data$eqtls$type <- 'quant'

# Perform colocalization
res <- coloc.abf(dataset1=data$gwas, dataset2=data$eqtls)

datasets <- list(datasets = list(meta_data = data$meta_data, gwas = data$gwas, eqtls = data$eqtls, coloc = res$results[,c("snp","SNP.PP.H4")], posterior = res$summary))

if (res$summary[6] >= 0.8){
  write(toJSON(datasets), paste("storage/out_files/", file_id, "/", gene, ".json", sep=''))
}
