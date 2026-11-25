// Annotate each patient's VCF with SnpEff.

process SNPEFF_ANN {
    tag "${patient} (${snpeff_genome})"
    publishDir "${params.outdir}/results/snpeff", mode: 'copy'

    input:
    tuple val(patient), path(vcf), val(snpeff_genome)

    output:
    tuple val(patient), path("${vcf.baseName}_ann.vcf"), emit: vcf
    tuple val(patient), path("${vcf.baseName}.fasta"),   emit: fasta

    script:
    def java_mem = task.memory ? Math.max(1, task.memory.toGiga() - 1) : 4
    """
    singularity exec --bind ${params.bind_paths} ${params.snpeff_sif} \\
        java -Xmx${java_mem}g -jar /snpEff/snpEff.jar ann \\
        -dataDir ${params.snpeff_datadir} \\
        -noStats \\
        -strict \\
        -hgvs1LetterAa \\
        -hgvs \\
        -canon \\
        -fastaProt ${vcf.baseName}.fasta \\
        ${snpeff_genome} \\
        ${vcf} \\
        > ${vcf.baseName}_ann.vcf

    # SnpEff writes no -fastaProt file when no variant changes a protein
    [ -e ${vcf.baseName}.fasta ] || touch ${vcf.baseName}.fasta
    """
}

workflow SNPEFF_ANNOTATION {
    take:
    ch_vcf   // tuple(patient, vcf, snpeff_genome)

    main:
    SNPEFF_ANN(ch_vcf)

    emit:
    vcf   = SNPEFF_ANN.out.vcf
    fasta = SNPEFF_ANN.out.fasta
}
