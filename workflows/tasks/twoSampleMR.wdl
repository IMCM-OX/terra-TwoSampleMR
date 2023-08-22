version 1.0

task twosamplemr {

    input {
        File exposure_summ_tab
        File outcome_summ_tab
        String marker
        String? npval = "5e-08"
        Float? clumping = 0.01
        String docker
    }

    parameter_meta {
        exposure_summ_tab: "Path to the 'exposure summary statistics' table. Tab-delimited file with eight required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval, Marker."
        outcome_summ_tab: "Path to the 'outcome summary statistics' table. Tab-delimited file with seven required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval."
        marker: "List of markers (comma seperated) to process through the Rscript"
        npval: "[Optional] nPval exposure filter. Default: '5e-08'"
        clumping: "[Optional] nClumping filter. Default: 0.01"
    }

    Int memory_mb = ceil(size(outcome_summ_tab, "MiB")) + 5000
    Int disk_size_gb = ceil(size(outcome_summ_tab, "GiB")) + 5

    command <<<
        set -euo pipefail
        Rscript /scripts/two_sample_mr__script_1.R \
            --exposures ~{marker} \
            --outstats ~{outcome_summ_tab} \
            --expstats ~{exposure_summ_tab} \
            --npval ~{npval} \
            --clumping ~{clumping}
    >>>

    output {
        File rds = "marker_set.mr_outputs_single_exposure.rds"
    }

    runtime {
        docker: "~{docker}"
        memory: "~{memory_mb} MiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        returnCodes: "*"
        continueOnReturnCode: true
    }
}
