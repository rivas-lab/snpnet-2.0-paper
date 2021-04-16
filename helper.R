computeStats <- function(pfile, ids, path, writefile=F) {
  keep_f       <- paste0(path, '.keep')
  gcount_tsv_f <- paste0(path, '.gcount.tsv')

    # To run plink2 --geno-counts, we write the list of IDs to a file
    data.frame(ID = ids) %>%
    tidyr::separate(ID, into=c('FID', 'IID'), sep='_') %>%
    data.table::fwrite(keep_f, sep='\t', col.names=F)

    # Run plink2 --geno-counts
    cmd_plink2 <- paste(
        'plink2',
        '--threads', 16,
        '--pfile', pfile,
        '--keep', keep_f,
        '--out', path,
        '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs'
    )

    system(cmd_plink2, intern=F, wait=T)

    # read the gcount file
    gcount_df <-
    data.table::fread(paste0(path, '.gcount')) %>%
    dplyr::mutate(
        stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
        stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
        stats_SDs   = stats_msts - stats_means * stats_means
    )

  out <- list()
  out[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
  out[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
  out[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)
  names(out[["means"]]) <- gcount_df %>% dplyr::pull(ID)
  to_exclude <- names(out[["means"]])[(out[["pnas"]] > 0.1) | (out[["means"]] < 2 * 0.001) | !(gcount_df[["#CHROM"]] %in% c(as.character(1:22), "X"))]

  if(writefile){
      data.frame(to_exclude) %>% data.table::fwrite(paste0(path, "var.exclude"), sep='\t', col.names=F)
  }
    return(list(keep=keep_f, exclude=paste0(path, "var.exclude"), df=gcount_df))
}