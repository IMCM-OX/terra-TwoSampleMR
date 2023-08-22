version 1.0

task generate_plot {

    input {
        File mr_sin_all_exposures
        File mr_sin_single_snp_exposures
        Float? fdr = 0.05
        Int? nsplit = 10
        String docker
    }

    parameter_meta {
        mr_sin_all_exposures: "Path to the mr_sin_all_exposures Rds object"
        mr_sin_single_snp_exposures: "Path to the mr_sin_single_snp_exposures Rds object"
        fdr: "[Optional] Threshold for the fdr adjusted p values used for the analyses. Default: 0.05"
        nsplit: "[Optional] The maximum number of exposures shown in each page of the pdf output file. Default: 10"
    }

    Int memory_mb = ceil(size(mr_sin_all_exposures, "MiB")) + 5000
    Int disk_size_gb = ceil(size(mr_sin_all_exposures, "GiB")) + 5
    
    command <<<
        set -euo pipefail
        Rscript /scripts/nv_wf1_multiANDsingle_v1.r \
            --msin ~{mr_sin_all_exposures} \
            --ssin ~{mr_sin_single_snp_exposures} \
            --fdr ~{fdr} \
            --nsplit ~{nsplit}
    >>>

    output {
        File multi_snp_joined_plot = "MR_scatter_forestplot_multi_snp_joined.pdf"
        File single_snp_joined_plot = "MR_scatter_forestplot_single_snp_joined.pdf"
        File multi_snp_table = "MR_scatter_forestplot_multi_snp_data.txt"
        File single_snp_table = "MR_scatter_forestplot_single_snp_data.txt"
    }
    
    runtime {
        docker: "~{docker}"
        memory: "~{memory_mb} MiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        returnCodes: "*"
        continueOnReturnCode: true
    }
}
