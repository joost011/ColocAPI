# packages
library(jsonlite)
library(coloc)
library(rjson)

args <- commandArgs(trailingOnly=TRUE)
file_path <- args[1]
ld_file_path <- args[2]
snp_list_file <- args[3]

# read files
data <- fromJSON(file = file_path)
LD <- read.table(file = ld_file_path)
snp_names <- read.table(file = snp_list_file)

# LD matrix
colnames(LD) <- unlist(snp_names)
rownames(LD) <- unlist(snp_names)

# attach GWAS data
data$gwas$type <- "cc"
data$gwas$N <- 276020
data$gwas$s <- 0.1284
data$gwas$LD <- as.matrix(LD)
# data$gwas$z <- data$eqtls$z

# attacth eQTL data
data$eqtls$type <- "quant"
data$eqtls$LD <- as.matrix(LD)
# data$eqtls$z <- data$eqtls$Zscore

# check datasets
# check_dataset(data$gwas, req=LD)
# check_dataset(data$eqtls, req=LD)
### Both returns NULL

# perform susie on GWAS and eQTL (get credible sets)
susie_gwas <- runsusie(data$gwas, n=data$gwas$N)
susie_eqtl <- runsusie(data$eqtls, n=31684)

summary(susie_gwas)
summary(susie_eqtl)
# res <- coloc.abf(dataset1=data$gwas, dataset2=data$eqtls)

# perform coloc on credible sets
susie.res <- coloc.susie(susie_gwas, susie_gwas)
susie.res$summary


if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.8",row=1,dataset1=data$gwas,dataset2=data$gwas)
  sensitivity(susie.res,"H4 > 0.8",row=2,dataset1=data$gwas,dataset2=data$gwas)
}
