workflow MevPuma {

    String identifier_choice

    # TODO: change these once we get miRNA priors
    Map[String, File] motifMap = {
        "Symbol":"s3://webmev-public/puma_mirna_priors.symbol.tsv",
        "Ensembl":"s3://webmev-public/puma_mirna_priors.ensg.tsv"
    }

    File motif_file = motifMap[identifier_choice]

    # A user uploaded exprs count matrix
    File exprs_file

    String output_name_prefix = "puma_matrix_results"

    call runPuma {
        input:
            motif_file = motif_file,
            exprs_file = exprs_file
    }

    output {
        File puma_output_matrix = runPuma.puma_output_matrix
    }
}

task runPuma {
    File motif_file
    File exprs_file

    # Don't want to expose this as user-defined variable, but want the 
    # flexibility to easily change it without a container rebuild.
    Int nmax = 25000

    String output_name = "puma_network.tsv"
    Int disk_size = 40

    command {
        python3 /opt/software/puma.py \
            --motif ${motif_file} \
            --output ${output_name} \
            --nmax ${nmax} \
            ${exprs_file}
    }

    output {
        File puma_output_matrix = "${output_name}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-puma"
        cpu: 16
        memory: "120 G"
        disks: "local-disk " + disk_size + " HDD"
    }
}
