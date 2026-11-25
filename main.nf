#!/usr/bin/env nextflow
// Mutation curation with bam-readcount.
// Samplesheet columns: patient,sample_name,status,bam_file,reference
// (reference is a genome version from params.genomes, e.g. hg38; blank = params.genome)

include { BAMREADCOUNT        } from './subworkflows/BAMREADCOUNT'
include { BAMRC_TO_VCF_WF     } from './subworkflows/BAMRC_TO_VCF'
include { SNPEFF_ANNOTATION   } from './subworkflows/SNPEFF_ANNOTATION'
include { SNPEFF_VCF_TO_TABLE } from './subworkflows/SNPEFF_VCF_TO_TABLE'

def findIndex(aln) {
    def p = aln.toString()
    def candidates = aln.name.endsWith('.cram')
        ? [ "${p}.crai", p.replaceFirst(/\.cram$/, '.crai') ]
        : [ "${p}.bai",  p.replaceFirst(/\.bam$/,  '.bai')  ]
    def idx = candidates.collect { c -> file(c) }.find { f -> f.exists() }
    if( !idx ) error "No index (.bai/.crai) found for ${aln}"
    return idx
}

// Resolve a samplesheet "reference" value (genome version, e.g. hg38 or GRCh38;
// blank = params.genome) to its params.genomes key, case-insensitively.
def genomeName(reference, sample) {
    def name = reference ?: params.genome
    def lookup = [:]
    params.genomes.keySet().each { k -> lookup[k.toString().toLowerCase()] = k.toString() }
    params.genome_aliases.each { alias, target -> lookup[alias.toString().toLowerCase()] = target.toString() }
    def key = lookup[name.toString().toLowerCase()]
    if( !key || !params.genomes[key] )
        error "Unknown reference '${name}' for ${sample}. Use a genome version from " +
              "${params.genomes.keySet() + params.genome_aliases.keySet()} or add it to params.genomes."
    return key
}

workflow {
    ch_samples = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .toList()
        .flatMap { rows -> (0..<rows.size()).collect { i -> rows[i] + [idx: i] } }
        .map { row ->
            def status = row.status?.toLowerCase()
            if( !(status in ['tumor', 'normal']) )
                error "Invalid status '${row.status}' for ${row.sample_name} (expected tumor/normal)"
            def bam   = file(row.bam_file, checkIfExists: true)
            def genome = genomeName(row.reference, row.sample_name)
            def fasta  = file(params.genomes[genome].fasta, checkIfExists: true)
            def meta   = [patient: row.patient, sample: row.sample_name, status: status,
                          idx: row.idx, genome: genome, reference: fasta.toString()]
            def fai   = file("${fasta}.fai", checkIfExists: true)
            tuple(meta, bam, findIndex(bam), fasta, fai)
        }

    BAMREADCOUNT(ch_samples, file(params.region, checkIfExists: true))
    BAMRC_TO_VCF_WF(BAMREADCOUNT.out.bamrc)

    // one SnpEff database per patient, from the samples' genome (--snpeff_genome overrides)
    ch_snpeff_genome = ch_samples
        .map { meta, bam, bai, fasta, fai ->
            tuple(meta.patient, params.snpeff_genome ?: params.genomes[meta.genome].snpeff) }
        .unique()

    SNPEFF_ANNOTATION(BAMRC_TO_VCF_WF.out.vcf.join(ch_snpeff_genome, failOnDuplicate: true))
    SNPEFF_VCF_TO_TABLE(SNPEFF_ANNOTATION.out.vcf)
}
