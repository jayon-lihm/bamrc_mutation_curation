// Flatten each patient's SnpEff-annotated VCF into a TSV table.

process SNPEFF_TO_TABLE {
    tag "${patient}"
    publishDir "${params.outdir}/results/protein_coding_table", mode: 'copy'

    input:
    tuple val(patient), path(ann_vcf)

    output:
    tuple val(patient), path("${ann_vcf.baseName}.protein_coding.tsv"), emit: table

    script:
    """
    ${params.python} ${projectDir}/src/snpeff_vcf_to_table.py \\
        --patient ${patient} \\
        -i ${ann_vcf} \\
        -o ${ann_vcf.baseName}.protein_coding.tsv
    """
}

workflow SNPEFF_VCF_TO_TABLE {
    take:
    ch_ann_vcf   // tuple(patient, annotated vcf)

    main:
    SNPEFF_TO_TABLE(ch_ann_vcf)

    emit:
    table = SNPEFF_TO_TABLE.out.table
}
