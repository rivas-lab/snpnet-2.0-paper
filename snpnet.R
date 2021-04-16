library(snpnet)

allargs = commandArgs(trailingOnly=TRUE)
ind = as.integer(allargs[1])


# ind = 1
phenotype.file = "WHERE THE PHENOTYPE AND THE COVARIATES ARE STORED"
genotype.pfile = "Plink2's pgen format file's prefix"
pheno_list = c("HC269", "HC382", "INI21001","INI50", "130696", "130700")
fam_list = c("binomial", "binomial", "gaussian", "gaussian", "cox", "cox")
phenotype <- pheno_list[ind]
family <- fam_list[ind]
results.dir = "Result directory"
dir.create(results.dir)

if(family != "cox"){
  covariates <- c("age", "sex", paste0("PC", 1:10))
  status = NULL
} else {
  covariates <- c( "sex", paste0("PC", 1:10))
  status = paste0("coxnet_status_f.", phenotype, ".0.0")
  phenotype = paste0("coxnet_y_f.", phenotype, ".0.0")
  phenotype.file = "Survival phenotype file"
}


configs <- list(
  results.dir = results.dir,
  nCores = 16,
  num.snps.batch = 2000,
  nlams.init = 20,
  nlams.delta = 10,
  save = TRUE,
  niter = 100,
  glmnet.thresh = 10^(-7),
  prevIter = 0,
  use.glmnetPlus = TRUE,
  early.stopping = TRUE,
  verbose = TRUE
)

start = Sys.time()
fit_snpnet <- snpnet(
  genotype.pfile = genotype.pfile,
  phenotype.file = phenotype.file,
  phenotype = phenotype,
  covariates = covariates,
  family = family,
  split.col = "split",
  status.col = status,
  mem = 32000,
  configs = configs
)
duration = Sys.time() - start
savefinal = list(duration=duration, fitresult= fit_snpnet)
save(savefinal, file=file.path(results.dir, "pre-result.RData"))

gc()
gcount.dir <- file.path(results.dir, "meta/snpnet.train.gcount")
prediction <- snpnet::predict_snpnet(
  fit = fit_snpnet,
  saved_path=results.dir, 
  new_genotype_file=genotype.pfile,
  new_phenotype_file=phenotype.file, 
  phenotype=phenotype,
  gcount_path=gcount.dir,
  status_col=status,
  covariate_names=covariates,
  family=family,
  idx=which.max(fit_snpnet$metric.val),
  split_col='split',
  split_name = c('test')
)

tosave = list(duration=duration, pred = prediction)
save(tosave, file=file.path(results.dir, "output.RData"))

