library(bigsnpr)
library(dplyr)
#########################
# Mostly adapted from https://privefl.github.io/bigstatsr/articles/penalized-regressions.html
# File paths are removed in this public version
########################
genotype.pfile = "Plink2's pgen format file's prefix"
results.dir = "Result directory"
ncores = 16

# # convert pgen format that we use to bed matrix
# cmd_plink2 <- paste(
#         'plink2',
#         '--pfile', genotype.pfile,
#         '--out', file.path(results.dir, "genotypes"),
#         '--make-bed'
#     )
# system(cmd_plink2, intern=F, wait=T)


# Get phenotype ready
allargs = commandArgs(trailingOnly=TRUE)
ind = as.integer(allargs[1])

phenotype.file <-  "WHERE THE PHENOTYPE AND THE COVARIATES ARE STORED"
pheno_list = c("HC269", "HC382", "INI21001","INI50")
fam_list = c("binomial", "binomial", "gaussian", "gaussian")
phenotype <- pheno_list[ind]
family <- fam_list[ind]

covariates <- c("age", "sex", paste0("PC", 1:10))
configs=list(zstdcat.path='zstdcat', zcat.path='zcat')


psam_id = snpnet::readIDsFromPsam(paste0(genotype.pfile, '.psam'))
phe_master = snpnet::readPheMaster(phenotype.file, psam_id, family, covariates, phenotype, NULL, 'split', configs)
if(family == "binomial"){
    phe_master[[phenotype]] = phe_master[[phenotype]] - 1
}

phe_test = phe_master %>% filter(split == "test")
phe_train = phe_master %>% filter(split %in% c("train", "val")) 

response = phe_train[[phenotype]]
covs = as.matrix(phe_train %>% select(all_of(covariates)))
ind.train = match(phe_train$ID, psam_id)



bedfile = "Location of the BED matrix"

tmpfile <- tempfile()

t1 = Sys.time()
snp_readBed2(bedfile, backingfile = tmpfile, ind.row=ind.train, ncores=ncores)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))
t2 = Sys.time()
attach_time = t2 - t1


G  <- obj.bigSNP$genotypes
G_impute = snp_fastImputeSimple(G, method="mean2", ncores = ncores)
t3 = Sys.time()
impute_time = t3 - t2


FUN = if(family == "binomial") big_spLogReg else big_spLinReg

fit_start = Sys.time()
fit <- FUN(X=G_impute, response, 
            covar.train=covs,
            pf.covar=rep(0, ncol(covs)),
            alphas=1,
            K=4,
            power_scale = 0,
            n.abort = 2, nlam.min = 30,
            ncores=16)  
fit_end  = Sys.time()

fit_time = fit_end - fit_start 



ind.test = match(phe_test$ID, psam_id)
tmpfile <- tempfile()

snp_readBed2(bedfile, backingfile = tmpfile, ind.row=ind.test, ncores=ncores)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))



G  <- obj.bigSNP$genotypes
# Still not quite right, but the mean imputation for training should not look at the test set.
G_impute = snp_fastImputeSimple(G, method="mean2", ncores = ncores)


pred = predict(
    fit,
    G_impute,
    covar.row = as.matrix(phe_test %>% select(all_of(covariates))),
    proba = FALSE,
    ncores = ncores
)


response = phe_test[[phenotype]]
if(family == "binomial"){
    pred.obj <- ROCR::prediction(pred, factor(response))
    auc.obj <- ROCR::performance(pred.obj, measure = 'auc')
    metric <- auc.obj@y.values[[1]]
} else {
    metric = 1 - sum( (response - pred)^2)/sum( (response - mean(response))^2)
}
fname = file.path(results.dir, paste0(phenotype, ".RData"))
save_list = list(attach_time=attach_time, impute_time=impute_time, fit_time=fit_time, metric=metric)
save(save_list, file=fname)
