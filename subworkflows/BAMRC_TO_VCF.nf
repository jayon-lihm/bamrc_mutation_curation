// Merge each patient's bam-readcount files into one multi-sample VCF.

process BAMRC_TO_VCF {
    tag "${patient}"
    publishDir "${params.outdir}/results/vcf", mode: 'copy'

    input:
    tuple val(patient), val(samples), path(bamrc_files, stageAs: 'bamrc/*'), path(fasta), path(fasta_fai)

    output:
    tuple val(patient), path("${patient}.bamrc.vcf"), emit: vcf

    script:
    """
    ${params.python} ${projectDir}/src/convert_bamrc_to_vcf.py \\
        --patient ${patient} \\
        --bamrc ${bamrc_files} \\
        --samples ${samples.join(' ')} \\
        --reference ${fasta} \\
        -o ${patient}.bamrc.vcf
    """
}

workflow BAMRC_TO_VCF_WF {
    take:
    ch_bamrc   // tuple(meta, bamrc) per sample; meta has patient, sample, idx, reference

    main:
    ch_patients = ch_bamrc
        .map { meta, bamrc -> tuple(meta.patient, meta, bamrc) }
        .groupTuple()
        .map { patient, metas, files ->
            // keep sample columns in samplesheet order
            def order = (0..<metas.size()).toList().sort { i -> metas[i].idx }
            def refs  = metas.collect { m -> m.reference }.unique()
            if( refs.size() > 1 )
                error "Patient ${patient} has samples with different references: ${refs}"
            def fasta = file(refs[0])
            tuple(patient,
                  order.collect { i -> metas[i].sample },
                  order.collect { i -> files[i] },
                  fasta,
                  file("${fasta}.fai"))
        }

    BAMRC_TO_VCF(ch_patients)

    emit:
    vcf = BAMRC_TO_VCF.out.vcf
}
