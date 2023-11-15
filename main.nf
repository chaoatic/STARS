#! usr/bin/env nextflow

// Concatenate raw read sequences 
process concatenateFastq {
    publishDir "${params.result_dir}/sequence", mode: "copy"

    input:
    path data_dir

    output:
    path "${data_dir}.fastq.gz"

    script:
    """
    cat ${data_dir}/* > ${data_dir}.fastq.gz
    """
}

// Visualization of raw read sequences
process rawSequenceVisualization {
    publishDir "${params.result_dir}/visualization/${concatenate_fastq.simpleName}", mode: "copy"

    input:
    path concatenate_fastq

    output:
    path "*"

    script:
    """
    NanoPlot --threads ${params.threads} --no_static --fastq ${concatenate_fastq}
    """
}

// Filter raw read sequences
process filterFastq {
    publishDir "${params.result_dir}/sequence", mode: "copy"

    input:
    path concatenate_fastq

    output:
    path "${concatenate_fastq.simpleName}_filter.fastq.gz"

    script:
    """
    gunzip -c ${concatenate_fastq} | \
    chopper \
        --minlength ${params.min_length} \
        --maxlength ${params.max_length} \
        --quality ${params.min_qscore} \
        --threads ${params.threads} |\
    gzip > ${concatenate_fastq.simpleName}_filter.fastq.gz
    """
}

// Visualization of filter read sequences
process filterSequenceVisualization {
    publishDir "${params.result_dir}/visualization/${filter_fastq.simpleName}", mode: "copy"

    input:
    path filter_fastq

    output:
    path "*"

    script:
    """
    NanoPlot --threads ${params.threads} --no_static --fastq ${filter_fastq}
    """
}

// Taxonomic classification using Emu
process taxonomicClassification {
    publishDir "${params.result_dir}/taxonomy/", mode: "copy"

    input:
    path filter_fastq

    output:
    path abundance
    path min_abundance
    path unclassified

    script:
    """
    emu abundance \
        --type map-ont \
        --min-abundance ${params.min_abundance} \
        --db ${params.db_dir} \
        --threads ${params.threads} \
        --keep-counts \
        --output-unclassified \
        ${filter_fastq}
    """
}

workflow {
    Channel
        .fromPath("${params.data_dir}/*", type: 'dir', checkIfExists: true)
        .set { input_ch }
    concatenateFastq(input_ch)
    rawSequenceVisualization(concatenateFastq.out)
    filterFastq(concatenateFastq.out)
    filterSequenceVisualization(filterFastq.out)
    taxonomicClassification(filterFastq.out)

    combineFeatureTable
    importArtifact
    taxonomicBarPlot
    diversityAnalysis
}