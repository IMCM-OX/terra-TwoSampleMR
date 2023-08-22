version 1.0

task get_all_exposures {
    
    input {
        File exposure_summ_tab
        Int? chunk_size = 50
    }

    parameter_meta {
        exposure_summ_tab: "Path to the 'exposure summary statistics' table. Tab-delimited file with eight required columns: SNP, beta, se, effect_allele, other_allele, eaf, pval, Marker."
        chunk_size: "[Optional] The parallel processing of markers will occur in subsets of user defined chunk sizes. Default: 50" 
    }

    Int memory_mb = ceil(size(exposure_summ_tab, "MiB")) + 5000
    Int disk_size_gb = ceil(size(exposure_summ_tab, "GiB")) + 5

    command <<<
        set -euo pipefail
        # Input TSV file path
        tsv_file=~{exposure_summ_tab}

        # Extract the column index for the column named "Marker"
        column_index=$(head -n 1 "$tsv_file" | tr '\t' '\n' | grep -n "Marker" | cut -d ":" -f 1)

        # Extract the IDs from the "Marker" column and create a list of unique IDs
        unique_ids=$(cut -f $column_index "$tsv_file" | tail -n +2 | sort | uniq)

        echo "$unique_ids" | xargs -L ~{chunk_size} echo | tr ' ' ','
    >>>
    
    output {
        Array[String] all_unique_exposures = read_lines(stdout())
    }

    runtime {
        docker: "debian:stable-20230502-slim@sha256:1529cbfd67815df9c001ed90a1d8fe2d91ef27fcaa5b87f549907202044465cb" 
        memory: "~{memory_mb} MiB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}
