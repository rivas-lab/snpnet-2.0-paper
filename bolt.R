library(snpnet) # Need a few helper function to do preprocessing
library(dplyr)
source("helper.R")


allargs = commandArgs(trailingOnly=TRUE)
ind = as.integer(allargs[1])
# ind = 1


phenotype.file = "WHERE THE PHENOTYPE AND THE COVARIATES ARE STORED"
genotype.pfile = "Plink2's pgen format file's prefix"
pheno_list = c("HC269", "HC382", "INI21001","INI50")
fam_list = c("binomial", "binomial", "gaussian", "gaussian")
phenotype <- pheno_list[ind]
family <- fam_list[ind]
results.dir = "Result directory"
unlink(results.dir,recursive=TRUE)
dir.create(results.dir)

covariates <- c("age", "sex", paste0("PC", 1:10))
configs=list(zstdcat.path='zstdcat', zcat.path='zcat')
configs[['gcount.full.prefix']] = file.path(results.dir, "gcount")
psam_id = readIDsFromPsam(paste0(genotype.pfile, '.psam'))
phe_master = readPheMaster(phenotype.file, psam_id, family, covariates, phenotype, NULL, 'split', configs)

if(family == "binomial"){
    phe_master[[phenotype]] = phe_master[[phenotype]] - 1
}

phe_test = phe_master %>% filter(split == "test")
phe_master = phe_master %>% filter(split %in% c("train", "val"))

# Set writefile=TRUE to save a keep file (samples that will be used) and an exlude file (variants that will be excluded)
keep_exclude = computeStats(genotype.pfile, phe_master$ID, configs[['gcount.full.prefix']], writefile=TRUE)


# Create BED files
cmd_plink2 <- paste(
        'plink2',
        '--pfile', genotype.pfile,
        '--keep', keep_exclude[["keep"]],
        '--exclude', keep_exclude[["exclude"]],
        '--out', file.path(results.dir, "genotypes"),
        '--make-bed'
    )
system(cmd_plink2, intern=F, wait=T)

## Genotype files writen at 
geno = file.path(results.dir, "genotypes")

##
pheno.file = file.path(results.dir, "phenotype")
data.table::fwrite(phe_master %>% select(all_of(c("FID", "IID", phenotype, covariates))), pheno.file, sep=' ', col.names=T)

##
cmd_bolt <- paste0(
    'bolt',
    ' --bfile=', geno,
    ' --phenoFile=', pheno.file,
    ' --phenoCol=', phenotype,
    ' --covarFile=', pheno.file,
    ' --qCovarCol=sex --qCovarCol=age --qCovarCol=PC{1:10}',
    ' --numThreads=16',
    ' --predBetasFile=', file.path(results.dir, "Betas"),
    ' --lmm',
    ' --geneticMapFile=WHERE_BOLT_IS_INSTALLED/BOLT-LMM_v2.3.5/tables/genetic_map_hg19_withX.txt.gz',
    ' --LDscoresFile=WHERE_BOLT_IS_INSTALLED/BOLT-LMM_v2.3.5/tables/LDSCORE.1000G_EUR.tab.gz',
    ' --statsFile=', file.path(results.dir, 'out')
)
start = Sys.time()
system(cmd_bolt, intern=F, wait=T)

end = Sys.time()
fit_time = end - start


beta_snp = data.table::fread(file.path(results.dir, "Betas"))



df = keep_exclude[[3]]
exclude_vec = as.data.frame(data.table::fread(keep_exclude[["exclude"]], header=F))[,1]
ind = which(! (df$ID %in% exclude_vec))

pgen_train = pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), sample_subset=match(phe_master$ID, psam_id))
Xtrain = pgenlibr::NewDense(pgen_train, ind) # Use my own imputation?
Xtrain_mean = pgenlibr::DenseTransMultv(Xtrain, rep(1/nrow(phe_master), nrow(phe_master)))

PRStrainmean = sum(Xtrain_mean * beta_snp$BETA)


pgen_test = pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), sample_subset=match(phe_test$ID, psam_id))
#Xtest = pgenlibr::NewDense(pgen_test, ind, df$stats_means[ind])
Xtest = pgenlibr::NewDense(pgen_test, ind, Xtrain_mean)

# Can't find a documentation, but the predBetas should be applied to centered variants?
PRStest = pgenlibr::DenseMultv(Xtest,beta_snp$BETA) - PRStrainmean


features_test = phe_test %>% select(all_of(covariates))
# The mean imputation should come from the training set

if(family == "gaussian"){
    # Get covariates coefficients
    glmmod <- stats::glm(
            stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates), collapse = " + "))),
            data = phe_master, family = family
        )

    alpha = coef(glmmod)[1]
    beta = c(coef(glmmod)[-1])
    pred = (alpha + as.matrix(features_test) %*% beta)[,1]
    response = phe_test[[phenotype]]

    pred_total = pred + PRStest
    metric = 1 - sum( (response - pred_total)^2)/sum( (response - mean(response))^2)
} else {
    # family is binomial
    PRStrain = pgenlibr::DenseMultv(Xtrain, beta_snp$BETA) - PRStrainmean
    phe_master$PRS = PRStrain
    glmwithPRS <- stats::glm(
            stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates, "PRS"), collapse = " + "))),
            data = phe_master, family = family
        )
    features_test = cbind(features_test, PRStest)
    alpha = coef(glmwithPRS)[1]
    beta = coef(glmwithPRS)[-1]
    pred = (alpha + as.matrix(features_test) %*% beta)[,1]

    pred.obj = ROCR::prediction(pred, factor(phe_test[[phenotype]]))
    auc.obj = ROCR::performance(pred.obj, measure = 'auc')
    metric = auc.obj@y.values[[1]]
}

fname = file.path(results.dir, paste0(phenotype, ".RData"))
save_list = list(fit_time=fit_time, metric=metric)
save(save_list, file=fname)


# Make sure that the ordering of beta matches that in the pgen files
if(!all(beta_snp$SNP == df$ID[ind])){
    stop("ordering not right")
}

if(!all(beta_snp$ALLELE1 == df$ALT[ind])){
    stop("encoding not right")
}
