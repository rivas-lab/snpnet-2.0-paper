# File paths are removed in this public version
library(snpnet)
library(dplyr)
library(survival)

allargs = commandArgs(trailingOnly=TRUE)
ind = as.integer(allargs[1])
# ind = 1

phenotype.file = "WHERE THE PHENOTYPE AND THE COVARIATES ARE STORED"
genotype.pfile = "Plink2's pgen format file's prefix"
pheno_list = c("HC269", "HC382", "INI21001","INI50", "130696", "130700")
fam_list = c("binomial", "binomial", "gaussian", "gaussian", "cox", "cox")
phenotype <- pheno_list[ind]
family <- fam_list[ind]


if(family != "cox"){
  covariates <- c("age", "sex", paste0("PC", 1:10))
  status = NULL
} else {
  covariates <- c( "sex", paste0("PC", 1:10))
  status = paste0("coxnet_status_f.", phenotype, ".0.0")
  phenotype = paste0("coxnet_y_f.", phenotype, ".0.0")
  phenotype.file = "COX PHENOTYPE FILE"
}



configs=list(zstdcat.path='zstdcat', zcat.path='zcat')
psam_id = readIDsFromPsam(paste0(genotype.pfile, '.psam'))
phe_master = readPheMaster(phenotype.file, psam_id, family, covariates, phenotype, status, 'split', configs)

phe_test = phe_master %>% filter(split == "test")
phe_train_val = phe_master %>% filter(split %in% c("train", "val")) 


if(family == "binomial"){
    phe_train_val[[phenotype]] = phe_train_val[[phenotype]] - 1
    phe_test[[phenotype]] = phe_test[[phenotype]] - 1
}
if(family != "cox"){
    glmmod <- stats::glm(
            stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates), collapse = " + "))),
            data = phe_train_val, family = family
        )
    response = phe_test[[phenotype]]
    pred_test = coef(glmmod)[1] + as.matrix(phe_test %>% select(all_of(covariates))) %*% coef(glmmod)[-1]
} else {
    form = as.formula(paste("Surv(", phenotype, ",", status, ") ~ ", paste(covariates, collapse = " + ")))
    coxmod = coxph(form, data=phe_train_val)
    pred_test =  as.matrix(phe_test %>% select(all_of(covariates))) %*% coef(coxmod)

}

if(family == "gaussian"){
     metric = 1 - sum( (response - pred_test)^2)/sum( (response - mean(response))^2)
} else if (family == "binomial"){
    pred.obj = ROCR::prediction(pred_test, factor(response - 1))
    auc.obj = ROCR::performance(pred.obj, measure = 'auc')
    metric = auc.obj@y.values[[1]]
} else {
    metric = cindex::CIndex(pred_test, phe_test[[phenotype]], phe_test[[status]])
}

print(phenotype)
print(metric)
