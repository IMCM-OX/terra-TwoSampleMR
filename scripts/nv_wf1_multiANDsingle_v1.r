#!/usr/bin/env Rscript

######################
## To be run as part of workflow one (aka the exporatory phase of a 2SMR analysis)
## Visualisation of the exposures with p adj values lower than a threshold
## through a single forest plot showing all methods and exposures

### Run this script by :
# Rscript nv_wf1_multiANDsingle_v1.r -m mr_sin_all_exposures.ID_test.rds  -s mr_sin_single_snp_exposures.ID_test_2.rds 
#######################/
# Main expected output: 
# 1. MR_scatter_forestplot_multi_snp_joined.pdf 
# 2. MR_scatter_forestplot_single_snp_joined.pdf
# 3. MR_scatter_forestplot_multi_snp_data.txt
# 4. MR_scatter_forestplot_single_snp_data.txt




## packages
suppressPackageStartupMessages(library("TwoSampleMR"))
library("ggplot2")
library("ggforestplot")
library("plyr")
library("optparse")
library("berryFunctions")
library("qpdf")

##
print("Checkmark 1: Installed and loaded packages")

## options
option_list = list(
    make_option(c("-m", "--msin"), default=NULL, help="Input of the mr_sin RDS file for all multi-SNP exposures",metavar="character"),
    make_option(c("-s", "--ssin"), default=NULL, help="Input of the mr_sin RDS file for all single-SNP exposures",metavar="character"),
    make_option(c("-f", "--fdr"), default=0.05, help="Threshold for the fdr adjusted p values used for the analyses",metavar="character"),
    make_option(c("-n", "--nsplit"), default=14, help="The maximum number of exposures shown in each page of the pdf output file",metavar="character")
    );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



if (is.null(opt$msin)){
    print_help(opt_parser)
    stop("Mandatory input missing: Input of the mr_sin RDS file for all multi-SNP exposures", call.=FALSE)
}
##
if (is.null(opt$ssin)){
    print_help(opt_parser)
    stop("Mandatory input missing: Input of the mr_sin RDS file for all single-SNP exposures", call.=FALSE)
}
##


fdrthreshold <- opt$fdr
n_forsplit<- opt$nsplit

############# Figure 1 for multi snp exposures #################

output_ms <- tryCatch(
    expr = {
        if (!(is.null(opt$msin))){
            print("Checkmark M1: Start of figure 1 for multi snp : preprocessing of data")

            msin_vis_wf1 <- readRDS(opt$msin)

            msin_vis_wf1 <- msin_vis_wf1[order(msin_vis_wf1$SNP),]

            msin_vis_wf1_padj<- msin_vis_wf1[which(msin_vis_wf1$SNP=="All - Inverse variance weighted"),] %>% 
            mutate(p_adj = p.adjust(p, "fdr"))
            rows_to_add <- (nrow(msin_vis_wf1)-nrow(msin_vis_wf1_padj))
            msin_vis_wf1_padj<- addRows(msin_vis_wf1_padj, n=rows_to_add)

            msin_vis_wf1_padj<- cbind(msin_vis_wf1, msin_vis_wf1_padj$p_adj)
            colnames(msin_vis_wf1_padj)[10]<- "p_adj"

            # Top 10 based on fdr adjusted p values:
            msin_vis_wf1_padj<- msin_vis_wf1_padj[order(msin_vis_wf1_padj$p_adj),]
            #print("Top 10 are:")
            #print(head(msin_vis_wf1_padj, 10))


            #fdrthreshold = 2.8  ############### !!!!! THIS SHOULD BE COMMENTED OUT  when running the non-test analysis (0.85 is the p adj value that gives 8 proteins)

            listpassingfdr_ms <- subset(msin_vis_wf1_padj, msin_vis_wf1_padj$p_adj< fdrthreshold) 
            listpassingfdr_ms <- unique(listpassingfdr_ms$exposure)

            #### If any of them pass fdr: #####
            if (length(listpassingfdr_ms)>0) {


                msin_vis_wf1_padj_f<- msin_vis_wf1_padj[which(msin_vis_wf1_padj$exposure%in%listpassingfdr_ms),]
                msin_vis_wf1_padj_f<- msin_vis_wf1_padj_f[order(msin_vis_wf1_padj_f$p_adj),]
               #


                print("Checkmark M2: creating the forest plot object")

                mr_forest_plots_sin <- mr_forest_plot(msin_vis_wf1_padj_f)

                for (i in 1:length(mr_forest_plots_sin)){
                    if (i==1) {
                    mr_forest_all_v2 <- subset(mr_forest_plots_sin[[i]][[1]], grepl("All", mr_forest_plots_sin[[i]][[1]]$SNP))
                    } else {
                    mr_forest_all_v2 <- rbind(mr_forest_all_v2, subset(mr_forest_plots_sin[[i]][[1]], grepl("All", mr_forest_plots_sin[[i]][[1]]$SNP)) ) 
                    }
                }
                
                for (i in 1:nrow(mr_forest_all_v2)) {
                    if (mr_forest_all_v2[i,6]!="All - Inverse variance weighted"){
                        mr_forest_all_v2[i, 10] <- c(".")
                    }
                }


                df_wf1_f1_ms <-
                    mr_forest_all_v2 %>%
                    mutate(
                    SNP = factor(
                        SNP,
                        levels =  c("All - Weighted mode","All - Weighted median","All - Simple mode","All - Maximum likelihood", "All - MR Egger", "All - Inverse variance weighted")
                    )
                    ) %>%
                    tidyr::drop_na(b)


                df_wf1_f1_ms<- arrange(df_wf1_f1_ms, desc(SNP), p_adj)
                print(head(df_wf1_f1_ms$p_adj, 20))
                print(head(df_wf1_f1_ms, 20))
                df_wf1_f1_ms$p_adj<- as.numeric(df_wf1_f1_ms$p_adj)
                df_wf1_f1_ms$p_adj<- format(df_wf1_f1_ms$p_adj, scientific=TRUE, nsmall=3, digits=2)
                df_wf1_f1_ms$p_adj[df_wf1_f1_ms$p_adj=="     NA"] <- c("     .")
                print(head(df_wf1_f1_ms$p_adj, 20))
                print(head(df_wf1_f1_ms, 20))
                
                


                print("Checkmark M3: creating the content of the pdf object")
                
                

                split_exposures_forest <- split(unique(df_wf1_f1_ms$exposure), ceiling(seq_along(unique(df_wf1_f1_ms$exposure))/n_forsplit))
                
                for (i in 1:length(split_exposures_forest)) {
                df_wf1_f1_ms_split<- df_wf1_f1_ms[which(df_wf1_f1_ms$exposure%in%split_exposures_forest[[i]]),]

                    print(head(df_wf1_f1_ms_split))

                    if (length(unique(df_wf1_f1_ms_split$exposure))>1){
                        y_plot= length(unique(df_wf1_f1_ms_split$exposure))+ 0.8
                        h_plot = length(unique(df_wf1_f1_ms_split$exposure))-1
                    } else {
                        y_plot = 10
                        h_plot = 3
                    }

                                        
                    # Draw a forestplot of betas for each split
                    fr_wadjp_i<- ggforestplot::forestplot(
                        df = df_wf1_f1_ms_split,
                        name = exposure,
                        estimate = b,
                        se = se,
                        colour = SNP,
                        pvalue = p,
                        xlab = "beta (95% CI)",
                        title =  paste0("Forest plot for multiple snp exposures (page: ", i, "/", length(split_exposures_forest), ")"),
                        logodds = FALSE
                    ) + coord_cartesian(xlim = c(-2, 2), clip="off" 
                    )+ geom_text(x = 2.02, 
                                hjust = 0,
                                size = 3, label= df_wf1_f1_ms_split$p_adj                
                    )+ annotate("text", x=2.15, y=y_plot, label= "FDR adjusted p value for IVW",               
                    ) + theme(plot.margin = unit(c(1,4,1,1), "lines"), legend.position = "right") # This widens the right margin

                    num_plot <- sprintf("%03d", as.numeric(i))


                    pdf(paste0("MR_scatter_forestplot_multi_snp_", num_plot, ".pdf"), height= h_plot, width = 15)
                    print(fr_wadjp_i)
                    dev.off()

                }
                
                
                print("Checkmark M4: End of figure 1 generation_when not emtpy for multi snps")

                qpdf::pdf_combine(input = list.files(pattern="MR_scatter_forestplot_multi_snp_*", full.names = TRUE), output = "MR_scatter_forestplot_multi_snp_joined.pdf")
                
                write.table(msin_vis_wf1_padj, file="MR_scatter_forestplot_multi_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)


                print("Checkmark M5: Combined the different pages of the figure 1 for multi snps")

            } else {  # and if none pass the fdr: 
            fr_wadjp<- c(0)
            h_plot=10

            
            write.table(fr_wadjp, file="MR_scatter_forestplot_multi_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)

            pdf("MR_scatter_forestplot_multi_snp_joined.pdf", height= h_plot, width = 15)
            print(fr_wadjp)
            dev.off()
            print("Checkmark M4: End of figure 1 generation_when emtpy for multi snps")

            }

            #
        }
    }, 
    error=function(e){
        fr_wadjp<- c(0)
        h_plot=10


        write.table(fr_wadjp, file="MR_scatter_forestplot_multi_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)

        pdf("MR_scatter_forestplot_multi_snp_joined.pdf", height= h_plot, width = 15)
        print(fr_wadjp)
        dev.off()
        print("Checkmark M2: End of figure 1 generation_when emtpy for single snp with error")

        message("Caught an error! Check if the input file for multiple snps is emtpy.")
        print(e)
    }
)

###



########## Figure 1 for single snp exposures #####################


output_ss <- tryCatch(
    expr = {
            if (!(is.null(opt$ssin))){
            print("Checkmark S1: Start of figure 1 for single snp : preprocessing of data")

            ssin_vis_wf1 <- readRDS(opt$ssin)

           
            ssin_vis_wf1 <- subset(ssin_vis_wf1, !(grepl("All", ssin_vis_wf1$SNP)))
            ssin_vis_wf1 <- subset(ssin_vis_wf1, !(grepl("No available data", ssin_vis_wf1$SNP)))
            ssin_vis_wf1_padj <- ssin_vis_wf1 %>% 
            mutate(p_adj = p.adjust(p, "fdr"))
            colnames(ssin_vis_wf1_padj)[10]<- "p_adj"

            # Top 10 based on fdr adjusted p values:
            ssin_vis_wf1_padj<- ssin_vis_wf1_padj[order(ssin_vis_wf1_padj$p_adj),]
            print("Top 10 are:")
            print(head(ssin_vis_wf1_padj, 10))

            #fdrthreshold= 2.8 ## 0.85 is the p adj value that gives 8 proteins

            listpassingfdr_ss<- subset(ssin_vis_wf1_padj, ssin_vis_wf1_padj$p_adj< fdrthreshold) 
            listpassingfdr_ss<- unique(listpassingfdr_ss$exposure)

            #print(listpassingfdr_ss)


            if (length(listpassingfdr_ss)>0) {
                ssin_vis_wf1_padj_f<- ssin_vis_wf1_padj[which(ssin_vis_wf1_padj$exposure%in%listpassingfdr_ss),]
                ssin_vis_wf1_padj_f<- ssin_vis_wf1_padj_f[order(ssin_vis_wf1_padj_f$p_adj),]
                ssin_vis_wf1_padj_f<- arrange(ssin_vis_wf1_padj_f, p_adj)
                #ssin_vis_wf1_padj_f$p_adj<- format(ssin_vis_wf1_padj_f$p_adj, scientific=TRUE, nsmall=3, digits=2)

                print("Checkmark S2: creating the forest plot object for single snp")
                
                df_wf1_f1_ss <-
                    ssin_vis_wf1_padj_f %>%
                    tidyr::drop_na(b)

                df_wf1_f1_ss <- arrange(df_wf1_f1_ss, p_adj)

                df_wf1_f1_ss$p_adj<- format(df_wf1_f1_ss$p_adj, scientific=TRUE, nsmall=3, digits=2)



                print("Checkmark S3: creating the content of the pdf object for single snp")
                
                split_exposures_forest <- split(unique(df_wf1_f1_ss$exposure), ceiling(seq_along(unique(df_wf1_f1_ss$exposure))/n_forsplit))
                
                for (i in 1:length(split_exposures_forest)) {
                    df_wf1_f1_ss_split<- df_wf1_f1_ss[which(df_wf1_f1_ss$exposure%in%split_exposures_forest[[i]]),]

                    print(head(df_wf1_f1_ss_split))

                    if (length(unique(df_wf1_f1_ss_split$exposure))>1){
                        y_plot= length(unique(df_wf1_f1_ss_split$exposure))+ 0.8
                        h_plot = length(unique(df_wf1_f1_ss_split$exposure))-1
                    } else {
                        y_plot = 10
                        h_plot = 3
                    }
                    # Draw a forestplot of betas for each split
                    fr_wadjp_i<- ggforestplot::forestplot(
                        df = df_wf1_f1_ss_split,
                        name = exposure,
                        estimate = b,
                        se = se,
                        pvalue = p,
                        xlab = "beta (95% CI)",
                        title =  paste0("Forest plot for multiple snp exposures (page: ", i, "/", length(split_exposures_forest), ")"),
                        logodds = FALSE
                    ) + coord_cartesian(xlim = c(-2, 2), clip="off"
                    )+ geom_text(x = 2.02,
                                hjust = 0,
                                size = 3, label= df_wf1_f1_ss_split$p_adj                
                    )+ annotate("text", x=1.8, y=y_plot, label= "FDR adjusted p value for Wald Ratio", 
                                
                    ) + theme(plot.margin = unit(c(1,4,1,1), "lines"), legend.position = "right") # This widens the right margin

                    num_plot <- sprintf("%03d", as.numeric(i))


                    pdf(paste0("MR_scatter_forestplot_single_snp_", num_plot, ".pdf"), height= h_plot, width = 15)
                    print(fr_wadjp_i)
                    dev.off()

                }

                print("Checkmark S4: End of figure 1 generation_when not emtpy for single snp")
                qpdf::pdf_combine(input = list.files(pattern="MR_scatter_forestplot_single_snp_*", full.names = TRUE), output = "MR_scatter_forestplot_single_snp_joined.pdf")
                
               
                write.table(ssin_vis_wf1_padj, file="MR_scatter_forestplot_single_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)

                print("Checkmark S5: Combining the dif pages of the figure 1 for single snp")  

            } else {
                fr_wadjp<- c(0)
                h_plot=10
                
                write.table(fr_wadjp, file="MR_scatter_forestplot_single_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)

                pdf("MR_scatter_forestplot_single_snp_joined.pdf", height= h_plot, width = 15)
                print(fr_wadjp)
                dev.off()
                print("Checkmark S4: End of figure 1 generation_when emtpy for single snp")

            }

        
        }
    }, 
    error= function(e){
        fr_wadjp<- c(0)
        h_plot=10

        
        write.table(fr_wadjp, file="MR_scatter_forestplot_single_snp_data.txt", quote= FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)

        pdf("MR_scatter_forestplot_single_snp_joined.pdf", height= h_plot, width = 15)
        print(fr_wadjp)
        dev.off()
        print("Checkmark S2: End of figure 1 generation_when emtpy for single snp with error")

        message("Caught an error! Check if the input file for single snps is emtpy.")
        print(e)
    }
)