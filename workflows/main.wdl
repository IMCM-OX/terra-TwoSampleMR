version 1.0

## This WDL pipeline executes 4 tasks:
## get_all_exposures : The first step in WDL workflow that extracts the list of exposures from the exposure summary statistics file provided as an input.
## twosamplemr : Second step in the WDL script that is executed using scatter function. It executes the Rscript mentioned below for every exposure ID obtained from Task 1 in parallel. The Rscript takes 3 inputs an exposure summary statistics file, an outcome summary statistics file, and a list of exposures to test.
## combine_objects : The third task in the workflow that collects the output files generated in the Task 2 and executes an Rscript that takes all the objects generated in the previous step and combines them into single .rds file.
## generate_plot : The final step that generates the forest plots.It calls a custom Rscript which is designed to generate two different plots, one for multi-SNPs and another for single-SNPs , and save them in the PDF format.
## The inputs are compatible with the Terra.bio platform. The details are provided under the `parameter_meta` sections.

import "./tasks/getAllExposures.wdl" as get_exp
import "./tasks/twoSampleMR.wdl" as tsmr
import "./tasks/combineObjects.wdl" as comb_obj
import "./tasks/generatePlot.wdl" as plt



workflow main {

    String pipeline_version = "1.0.0"
    String container_src = "ghcr.io/imcm-ox/twosamplemr:~{pipeline_version}"

    input {
        File exposure_summary_statistics
        File outcome_summary_statistics
    }

    parameter_meta {
        exposure_summary_statistics: "Path to the 'exposure summary statistics' table. (Tab-delimited file with eight required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval, Marker)."
        outcome_summary_statistics: "Path to the 'outcome summary statistics' table. (Tab-delimited file with seven required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval)."
    }

    call get_exp.get_all_exposures {
        input: exposure_summ_tab = exposure_summary_statistics
    }
    
    scatter (id_set in get_all_exposures.all_unique_exposures) {
        call tsmr.twosamplemr {
            input: exposure_summ_tab = exposure_summary_statistics, outcome_summ_tab = outcome_summary_statistics, marker = id_set, docker = container_src
        }
    }

    #Returns: an Array of all non-None values in the input array.
    Array[File] rds_objects = select_all(twosamplemr.rds)

    call comb_obj.combine_objects {
        input: rds_objs = rds_objects, docker = container_src
    }

    call plt.generate_plot {
        input: mr_sin_all_exposures = combine_objects.mr_sin_all_exposures, mr_sin_single_snp_exposures = combine_objects.mr_sin_single_snp_exposures, docker = container_src
    }

    output {
        File multi_snp_joined_plot = generate_plot.multi_snp_joined_plot
        File single_snp_joined_plot = generate_plot.single_snp_joined_plot
        File multi_snp_table = generate_plot.multi_snp_table
        File single_snp_table = generate_plot.single_snp_table
        File mr_sin_all_exposures = combine_objects.mr_sin_all_exposures
        File mr_sin_single_snp_exposures = combine_objects.mr_sin_single_snp_exposures
        File mr_res_all_exposures = combine_objects.mr_res_all_exposures
        File dat_all_exposures = combine_objects.dat_all_exposures
        File mr_loo_all_exposures = combine_objects.mr_loo_all_exposures
    }
    
    meta {
        description: "A WDL-based workflow that utilizes the R package 'TwoSampleMR' to generate single SNP and multi SNP forest plots. This workflow is designed to facilitate Mendelian randomization analysis using TwoSampleMR package and visualize the results in the form of forest plots."
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }
}
