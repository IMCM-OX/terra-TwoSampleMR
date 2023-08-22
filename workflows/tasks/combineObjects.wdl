version 1.0

task combine_objects {

    input {
        Array[File] rds_objs
        String docker
    }

    parameter_meta {
        rds_objs: "Files of file paths to the rds objects (comma seperated) generated in the previous step"
    }

    Int memory_mb = ceil(size(rds_objs, "MiB")) + 5000
    Int disk_size_gb = ceil(size(rds_objs, "GiB")) + 5
    
    command <<<
        set -euo pipefail
        ulimit -u unlimited
        Rscript /scripts/gather_mr_outputs__script_2.R \
            --mrresults ~{sep="," rds_objs} \
            --analysisid "combined"
    >>>

    output {
        File mr_sin_all_exposures = "mr_sin_all_exposures.combined.rds"
        File mr_sin_single_snp_exposures = "mr_sin_single_snp_exposures.combined.rds"
        File mr_res_all_exposures = "mr_res_all_exposures.combined.rds"
        File dat_all_exposures = "dat_all_exposures.combined.rds"
        File mr_loo_all_exposures = "mr_loo_all_exposures.combined.rds"
    }
    
    runtime {
        docker: "~{docker}"
        memory: "~{memory_mb} MiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        returnCodes: "*"
        continueOnReturnCode: true
    }
}
