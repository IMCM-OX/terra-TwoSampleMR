#!/usr/bin/env Rscript
# Libraries ------------------------------------------------------------------------------ #
library("optparse")
suppressPackageStartupMessages(library("TwoSampleMR"))

option_list = list(
    make_option(c("-p", "--expstats"), type="character", default=NULL, help="Input summary data for exposures", metavar="character"),
    make_option(c("-t", "--outstats"), type="character", default=NULL, help="Input summary statistics for outcome", metavar="character"), 
    make_option(c("-e", "--exposures"), type="character", default=NULL, help="Name of the Exposure to test. (It can also be a list of comma seperated exposures)", metavar="character"),
    make_option(c("-v", "--npval"), type="double", default=5e-08, help="nPval exposure filter. [default= %default]", metavar="character"),
    make_option(c("-c", "--clumping"), type="double", default=0.01, help="nClumping filter [default= %default]", metavar="character"), 
    make_option(c("-s", "--stub"), type="character", default="mr_outputs_single_exposure", help="Stub for output RDS files. [default= %default]", metavar="character") 
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$expstats)){
    print_help(opt_parser)
    stop("Input summary data for exposures is mandatory", call.=FALSE)
}

if (is.null(opt$outstats)){
    print_help(opt_parser)
    stop("Input summary statistics for outcome is mandatory", call.=FALSE)
}

if (is.null(opt$exposures)){
    print_help(opt_parser)
    stop("List of exposures to test is mandatory", call.=FALSE)
}


# PARAMETERS ------------------------------------------------------------------------- #
nPval_exposure = NULL
nClumping = NULL

if (!is.null(opt$npval)){
    nPval_exposure = opt$npval
}

if (!is.null(opt$clumping)){
    nClumping = opt$clumping
}
# ------------------------------------------------------------------------------------ #



###############################################################
#### MR analysis
###############################################################
# Upload the outcome summary statistics data with column names: SNP effect_allele other_allele beta se pval eaf
outcome_stats<-TwoSampleMR::read_outcome_data(
    filename=file.path(opt$outstats), sep="\t",
    snp_col = "SNP", beta_col = "beta", se_col="se", 
    effect_allele_col = "effect_allele", other_allele_col="other_allele",
    eaf_col="eaf", pval_col="pval"
)

# Upload the exposure summary statistics data save before with column names SNP, beta, se, effect_allele, other_allele, pval, eaf, Marker
exposure_stats<-TwoSampleMR::read_exposure_data(
    filename=file.path(opt$expstats), 
        sep="\t", 
        snp_col = "SNP", beta_col = "beta", se_col="se", 
        effect_allele_col = "effect_allele", other_allele_col = "other_allele",
        pval_col="pval", eaf_col="eaf", phenotype_col="Marker"
        )

exposure_stats_sig<-exposure_stats[which(exposure_stats$pval.exposure < nPval_exposure),]

exposures_array <- unique(strsplit(opt$exposures, ",")[[1]])


mr_outputs_lists <- list()
for( exposure in exposures_array ){

  tryCatch(
    {
      error_code <- -1  
      #exposure_stats_keep<-exposure_stats_sig[which(exposure_stats_sig$exposure%in%exposures_array),]

      exposure_stats_keep<-exposure_stats_sig[which(exposure_stats_sig$exposure == exposure),]

      # LD prunning
      exposure_stats_keep_pruned<-TwoSampleMR::clump_data(exposure_stats_keep, clump_r2=nClumping)

      # Extracting matching SNPs
      outcome_snps_in_exposure_snps <- outcome_stats$SNP %in% exposure_stats_keep_pruned$SNP
      outcome_stats_final_subset <- outcome_stats[which(outcome_snps_in_exposure_snps==TRUE),]
      num_snps <- sum(outcome_snps_in_exposure_snps)


      #suggest adding:
      # exposure_snps_in_outcome_snps <-exposure_stats_keep_pruned$SNP %in% outcome_stats$SNP
      # exposure_stats_final_subset <- exposure_stats_keep_pruned[which(exposure_snps_in_outcome_snps==TRUE),]

      # Harmonization
      dat<-TwoSampleMR::harmonise_data(exposure_dat= exposure_stats_keep_pruned, outcome_dat = outcome_stats_final_subset, 2)
      # suggest:
      #dat<-TwoSampleMR::harmonise_data(exposure_dat= exposure_stats_final_subset, outcome_dat = outcome_stats_final_subset, 2)
  

      mr_res <- TwoSampleMR::mr(dat)
      mr_sin <- TwoSampleMR::mr_singlesnp(dat)

      mr_het <- TwoSampleMR::mr_heterogeneity(dat)
      mr_plt <- TwoSampleMR::mr_pleiotropy_test(dat)
      mr_loo <- TwoSampleMR::mr_leaveoneout(dat, method= mr_ivw)
      
      error_code <- 0
    },
    error=function(e){
      message('An error occurred')
      print(e)
    },
    finally={
      if(error_code == 0){
        mr_outputs_list <- list(exposure,
			   num_snps,
			   dat,
			   mr_res,
			   mr_sin,
			   mr_het,
			   mr_plt,
			   mr_loo)
      }else{
        mr_outputs_list <- list(exposure,
			   -1,
			   data.frame(),
			   data.frame(),
			   data.frame(),
			   data.frame(),
			   data.frame(),
			   data.frame())
      }
      names(mr_outputs_list) <- c("exposure.id","nsnps","dat","mr_res","mr_sin","mr_het","mr_plt","mr_loo")

      mr_outputs_lists<-append(mr_outputs_lists,list(mr_outputs_list))

      #saveRDS(mr_outputs_list, file=file.path(getwd(),paste(exposure,opt$stub,"rds",sep=".")))
      
    }
  )
}
saveRDS(mr_outputs_lists, file=file.path(getwd(),paste("marker_set",opt$stub,"rds",sep=".")))
