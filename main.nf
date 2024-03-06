params.motif_files_map = [
    symbol: "s3://webmev-public/puma_mirna_priors.symbol.tsv",
    ensembl: "s3://webmev-public/puma_mirna_priors.ensg.tsv"
]
params.nmax = 25000

process run_puma {

    tag "Run PUMA"
    publishDir "${params.output_dir}/MevPuma.puma_output_matrix", mode:"copy", pattern:"${output_fname}"
    container "ghcr.io/web-mev/mev-puma"
    cpus 16
    memory '120 GB'

    input:
        path exprs_file
        path motif_file

    output:
        path "${output_fname}"

    script:
        output_fname = "puma_network.tsv"
        """
        python3 /usr/local/bin/puma.py \
            --motif ${motif_file} \
            --output ${output_fname} \
            --nmax ${params.nmax} \
            ${exprs_file}
        """
}

workflow {
    motif_file_ch = Channel.fromPath(params.motif_files_map[params.identifier_choice.toLowerCase()])
    run_puma(params.exprs_file, motif_file_ch)
}