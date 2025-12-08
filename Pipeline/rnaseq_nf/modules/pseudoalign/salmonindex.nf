// salmon index function

nextflow.enable.dsl=2

process SALMON_INDEX {

    tag "${transcriptome.simpleName}"

    publishDir "${params.outdir}/salmon_index", mode: 'copy', overwrite: true

    input:
        path transcriptome

    output:
        path "salmon_index"
    
    script:
    """
    salmon index -t ${transcriptome} -i salmon_index
    """

}
