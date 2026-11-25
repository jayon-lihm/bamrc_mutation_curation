// Run bam-readcount (SNVs + indels) over a region list for each BAM/CRAM.

process BAM_READCOUNT {
    tag "${meta.sample}"
    publishDir path: { "${params.outdir}/results/bamreadcount/${meta.patient}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bam_index), path(fasta), path(fasta_fai)
    path region

    output:
    tuple val(meta), path("${meta.sample}.bamrc.tsv"), emit: bamrc

    script:
    def idx_ext = bam.name.endsWith('.cram') ? 'crai' : 'bai'
    """
    # bam-readcount expects the index next to the alignment as <bam>.<ext>
    [ -e ${bam}.${idx_ext} ] || ln -s ${bam_index} ${bam}.${idx_ext}

    singularity exec --bind ${params.bind_paths} ${params.bamrc_sif} \\
        bam-readcount \\
        -f ${fasta} \\
        -l ${region} \\
        -i \\
        -b ${params.min_base_qual} \\
        -w ${params.max_warnings} \\
        ${bam} \\
        > ${meta.sample}.bamrc.tsv
    """
}

workflow BAMREADCOUNT {
    take:
    ch_samples   // tuple(meta, bam, bam_index, fasta, fasta_fai)
    region       // path to BED region list

    main:
    BAM_READCOUNT(ch_samples, region)

    emit:
    bamrc = BAM_READCOUNT.out.bamrc
}
