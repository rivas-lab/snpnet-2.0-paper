library(snpnet) # Need a few helper function to do preprocessing
library(dplyr)
library(data.table)
library(bigsnpr)

###############
# Heavily adapted from https://github.com/privefl/paper-ldpred2/blob/master/code/run-ldpred2-gwide.R#L34-L118
###############
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
#unlink(results.dir,recursive=TRUE)
dir.create(results.dir)

covariates <- c("age", "sex", paste0("PC", 1:10))
covariatesNoSex <- c("age", paste0("PC", 1:10))
configs=list(zstdcat.path='zstdcat', zcat.path='zcat')
configs[['gcount.full.prefix']] = file.path(results.dir, "gcount")
psam_id = readIDsFromPsam(paste0(genotype.pfile, '.psam'))
phe_master = readPheMaster(phenotype.file, psam_id, family, covariates, phenotype, NULL, 'split', configs) %>% filter(split %in% c('train', 'test', 'val'))

phe_test = phe_master %>% filter(split == "test")
phe_val = phe_master %>% filter(split == "val") 
phe_train = phe_master %>% filter(split == "train") 


phe_train %>% select(all_of(c("FID", "IID", phenotype, covariatesNoSex))) %>% data.table::fwrite(file=file.path(results.dir, "covs"), sep='\t')

cmd_plink2 <- paste(
        'plink2',
        '--pfile', genotype.pfile,
        '--chr 1-22',
        '--pheno', file.path(results.dir, "covs"),
        '--pheno-col-nums 3',
        '--covar',  file.path(results.dir, "covs"),
        '--covar-col-nums', paste0('4-', 3 + length(covariatesNoSex)), 
        '--glm sex',
        '--keep', file.path(results.dir, "covs"),
        '--covar-variance-standardize',
        '--out', file.path(results.dir, 'GWAS'),
        '--threads 16'
    )
t0 = Sys.time()
system(cmd_plink2, intern=F, wait=T)
Gwas_time = Sys.time() - t0
# Only need the effect of the snps
result_file = file.path(results.dir, list.files(results.dir, pattern="glm"))
if(family == "binomial"){
    ncontrol = sum(phe_train[[phenotype]] == 1)
    ncase = nrow(phe_train) - ncontrol
    sumstats = data.table::fread(result_file )%>% filter(TEST == "ADD" & (!is.na(OR)))  %>% select(all_of(c("#CHROM", "ID", "POS", "REF", "ALT", "A1", "OR", "LOG(OR)_SE", "P")))
    # REF will become A0
    tmp.ind = which(sumstats$A1 != sumstats$ALT) # in this case REF is a1, so alt is a0
    sumstats$REF[tmp.ind] = sumstats$ALT[tmp.ind]
    sumstats = sumstats %>% select(all_of(c("#CHROM", "ID", "POS", "REF", "A1", "OR", "LOG(OR)_SE", "P")))
    sumstats$OR = log(sumstats$OR)
    sumstats$n_eff = 4 / ((1 / ncontrol) + (1 / ncase))
    names(sumstats) =  c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
} else {
    sumstats = data.table::fread(result_file) %>% filter(TEST == "ADD" & (!is.na(BETA)))  %>% select(all_of(c("#CHROM", "ID", "POS", "REF", "ALT", "A1", "BETA", "SE", "P")))
    tmp.ind = which(sumstats$A1 != sumstats$ALT) # in this case REF is a1, so alt is a0
    sumstats$REF[tmp.ind] = sumstats$ALT[tmp.ind]
    sumstats = sumstats %>% select(all_of(c("#CHROM", "ID", "POS", "REF", "A1", "BETA", "SE", "P")))
    sumstats$n_eff = nrow(phe_train)
    names(sumstats) =  c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
}


ind.val = match(phe_val$ID, psam_id)
bimfile = read.table("BIM_FILE_FOLDER/genotypes.bim", header=FALSE)
bimfile = bimfile$V1
ind.var = which(bimfile %in% c(1:22))



bedfile = "BED_FILE_FOLDER/genotypes.bed"

tmpfile <- tempfile()
ncores = 16
t1 = Sys.time()
snp_readBed2(bedfile, backingfile = tmpfile, ind.row=ind.val, ind.col=ind.var, ncores=ncores)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))
t2 = Sys.time()
attach_time = t2 - t1


G  <- obj.bigSNP$genotypes
G_impute = snp_fastImputeSimple(G, method="mean2", ncores = ncores)
t3 = Sys.time()
impute_time = t3 - t2

obj.bigSNP$map$chromosome = as.integer(obj.bigSNP$map$chromosome)
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map)

tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)


t4 = Sys.time()
POS2 <- snp_asGeneticPos(CHR, POS, dir = "1000-genomes-genetic-maps/interpolated_OMNI", ncores = ncores)
for(chr in 1:22){
    print(chr)
    ind.chr <- which(info_snp$chr == chr)         
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

    corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = ncores, 
            infos.pos = POS2[ind.chr2], size = 3 / 1000)

    if(chr == 1){
        df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
        corr <- as_SFBM(corr0, tmp)
        ld <- Matrix::colSums(corr0^2)
    } else {
        df_beta = rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
        corr$add_columns(corr0, nrow(corr))
        ld <- c(ld, Matrix::colSums(corr0^2))
    }
}
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
t5 = Sys.time()
corr_time = t5 - t4


multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.9, 30),
                                 ncores = ncores)
t6 = Sys.time()
fit_time = t6 - t5
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G_impute, beta_auto,
                           ind.col = df_beta[["_NUM_ID_"]])
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])


# Evaluate beta


tmpfile2 <- tempfile()

snp_readBed2(bedfile, backingfile = tmpfile2, ind.row=match(phe_master$ID, psam_id), ind.col=ind.var, ncores=ncores)
# Attach the "bigSNP" object in R session
obj.bigSNP.all <- snp_attach(paste0(tmpfile2, ".rds"))



G.all  <- obj.bigSNP.all$genotypes

G_impute_all = snp_fastImputeSimple(G.all , method="mean2", ncores = ncores)

phe_train_val = phe_master %>% filter(split %in% c('train', 'val'))

pred_all = big_prodVec(G_impute_all, final_beta_auto, ind.col = df_beta[["_NUM_ID_"]], ncores=ncores)

pred_train_val = pred_all[match(phe_train_val$ID, phe_master$ID)]
pred_test = pred_all[match(phe_test$ID, phe_master$ID)]
# pred_test = big_prodVec(G_impute_all, final_beta_auto, ind.col = df_beta[["_NUM_ID_"]], ind.row=match(phe_test$ID, phe_master$ID), ncores=ncores)

phe_train_val$PRS = pred_train_val
phe_test$PRS = pred_test
if(family == "binomial"){
    phe_train_val[[phenotype]] = phe_train_val[[phenotype]] - 1
    glmmod = stats::glm(
            stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates, "PRS"), collapse = " + "))),
            data = phe_train_val, family = family
        )
    pred_final = coef(glmmod)[1] + as.matrix(phe_test %>% select(all_of(c(covariates, "PRS")))) %*% coef(glmmod)[-1]

    pred.obj <- ROCR::prediction(pred_final, factor(phe_test[[phenotype]]))
    auc.obj <- ROCR::performance(pred.obj, measure = 'auc')
    metric <- auc.obj@y.values[[1]]
} else {
     glmmod = stats::glm(
            stats::as.formula(paste(phenotype, " ~ ", paste(c(1, covariates, "PRS"), collapse = " + "))),
            data = phe_train_val, family = family
        )
    pred_final = coef(glmmod)[1] + as.matrix(phe_test %>% select(all_of(c(covariates, "PRS")))) %*% coef(glmmod)[-1]


    response = phe_test[[phenotype]]
    metric = 1 - sum( (response - pred_final)^2)/sum( (response - mean(response))^2)
}


fname = file.path(results.dir, paste0(phenotype, ".RData"))
save_list = list(Gwas_time=Gwas_time, attach_time=attach_time, impute_time=impute_time, fit_time=fit_time, metric=metric)
save(save_list, file=fname)

