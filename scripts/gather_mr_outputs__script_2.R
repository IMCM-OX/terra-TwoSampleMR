#!/usr/bin/env Rscript
# Libraries ------------------------------------------------------------------------------ #
library("optparse")
suppressPackageStartupMessages(library("dplyr"))

option_list = list(
    make_option(c("-r", "--mrresults"), type="character", default=NULL, help="Comma separated list of RDS filenames.", metavar="character"),
    make_option(c("-a", "--analysisid"), type="character", default=NULL, help="Analysis ID foroutput RDS files.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$mrresults)){
    print_help(opt_parser)
    stop("MR results for gathering are mandatory.", call.=FALSE)
}


if (is.null(opt$analysisid)){
    print_help(opt_parser)
    stop("Analysis ID for output RDS files is mandatory.", call.=FALSE)
}


# Read in the mr outputs lists for all exposures ---> These have been input as a comma separated list of RDS filenames *** 

all_mr_outputs_lists_aux <- lapply(unique(strsplit(opt$mrresults, ",")[[1]]), readRDS)

all_mr_outputs_lists<-unlist(all_mr_outputs_lists_aux, recursive=FALSE)


if (length(all_mr_outputs_lists)==0){
    stop("Error: List of MR results for gathering is empty. Please run script again, but with a populated list.", call.=FALSE)
}

# ***********************************************************************************************************************************

# Now make boolean indicator arrays to identify exposures with one (single), and more than one (multi) SNP. *************************
single_snp_exposures <- unlist(lapply(all_mr_outputs_lists,function(x) x$nsnps==1))
multi_snp_exposures <- unlist(lapply(all_mr_outputs_lists,function(x) x$nsnps>1))

# ***********************************************************************************************************************************

# ***********************************************************************************************************************************
# Save MR output data frames for single SNP exposures.
# ***********************************************************************************************************************************
dat_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$dat))[which(single_snp_exposures)])
mr_res_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_res))[which(single_snp_exposures)])
mr_sin_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_sin))[which(single_snp_exposures)])
# This is a dataframe where each exposure has a row, but results are all NA... Keeping it as possible required for figure generation.
mr_loo_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_loo))[which(single_snp_exposures)])

# These data frames are empty and are not used downstream in figure generation, therefore not making them.
#mr_het_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_het))[which(single_snp_exposures)])
#mr_plt_single_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_plt))[which(single_snp_exposures)])

saveRDS(dat_single_snp_exposures, file=file.path(getwd(),paste("dat_single_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_res_single_snp_exposures, file=file.path(getwd(),paste("mr_res_single_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_sin_single_snp_exposures, file=file.path(getwd(),paste("mr_sin_single_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_loo_single_snp_exposures, file=file.path(getwd(),paste("mr_loo_single_snp_exposures",opt$analysisid,"rds",sep=".")))

# ***********************************************************************************************************************************
# Save MR output data frames for multi SNP exposures.
# ***********************************************************************************************************************************
dat_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$dat))[which(multi_snp_exposures)])
mr_res_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_res))[which(multi_snp_exposures)])
mr_sin_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_sin))[which(multi_snp_exposures)])
mr_het_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_het))[which(multi_snp_exposures)])
mr_plt_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_plt))[which(multi_snp_exposures)])
mr_loo_multi_snp_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_loo))[which(multi_snp_exposures)])


saveRDS(dat_multi_snp_exposures, file=file.path(getwd(),paste("dat_multi_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_res_multi_snp_exposures, file=file.path(getwd(),paste("mr_res_multi_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_sin_multi_snp_exposures, file=file.path(getwd(),paste("mr_sin_multi_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_het_multi_snp_exposures, file=file.path(getwd(),paste("mr_het_multi_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_plt_multi_snp_exposures, file=file.path(getwd(),paste("mr_plt_multi_snp_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_loo_multi_snp_exposures, file=file.path(getwd(),paste("mr_loo_multi_snp_exposures",opt$analysisid,"rds",sep=".")))


# ***********************************************************************************************************************************
# Finally, in case they are neeeded, save MR output data frames for all exposures (regardless of number of SNPs) together.
# ***********************************************************************************************************************************
dat_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$dat)))
mr_res_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_res)))
mr_sin_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_sin)))
mr_het_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_het)))
mr_plt_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_plt)))
mr_loo_all_exposures <- bind_rows((lapply(all_mr_outputs_lists,function(x) x$mr_loo)))


saveRDS(dat_all_exposures, file=file.path(getwd(),paste("dat_all_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_res_all_exposures, file=file.path(getwd(),paste("mr_res_all_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_sin_all_exposures, file=file.path(getwd(),paste("mr_sin_all_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_het_all_exposures, file=file.path(getwd(),paste("mr_het_all_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_plt_all_exposures, file=file.path(getwd(),paste("mr_plt_all_exposures",opt$analysisid,"rds",sep=".")))
saveRDS(mr_loo_all_exposures, file=file.path(getwd(),paste("mr_loo_all_exposures",opt$analysisid,"rds",sep=".")))

