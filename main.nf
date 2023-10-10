#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    INPUT FILES
======================================================================================== */
params.input = null
input_files  = file(params.input)


/* ========================================================================================
    OUTPUT DIRECTORY
======================================================================================== */
params.outdir = false
if(params.outdir){
    outdir = params.outdir
} else {
    outdir = '.'
}

/* ========================================================================================
    SKIP STEPS
======================================================================================== */
params.skip_qc           = false
params.skip_fastq_screen = false


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.fastq_screen_conf = '/cluster/work/nme/software/config/fastq_screen.conf' // FastQ Screen config file directory
params.genome            = 'GRCm39'    // Default genome
params.strandness        = 'smartseq2' // Library strandedness
params.single_end        = false       // Force to input files to be single-end


/* ========================================================================================
    EXTRA PARAMETERS
======================================================================================== */
params.fastqc_args         = ''
params.fastq_screen_args   = ''
params.trim_galore_args    = ''
params.star_align_args     = ''
params.hisat2_args         = ''
params.featurecounts_args  = ''
params.multiqc_args        = ''
params.samtools_view_args  = ''
params.samtools_sort_args  = ''
params.samtools_index_args = ''

fastqc_args         = params.fastqc_args       
fastq_screen_args   = params.fastq_screen_args
trim_galore_args    = params.trim_galore_args 
star_align_args     = params.star_align_args 
hisat2_args         = params.hisat2_args 
featurecounts_args  = params.featurecounts_args
multiqc_args        = params.multiqc_args
samtools_view_args  = params.samtools_view_args
samtools_sort_args  = params.samtools_sort_args
samtools_index_args = params.samtools_index_args

/* ========================================================================================
    ALIGNER
======================================================================================== */
params.aligner = 'star'


/* ========================================================================================
    HISAT2 PARAMETERS
======================================================================================== */
// no_softclip
params.no_softclip = true

if(params.no_softclip){
    hisat2_args += " --no-softclip "
}

/* ========================================================================================
    STAR PARAMETERS
======================================================================================== */
// outSAMtype
params.outSAMtype_file = 'BAM'
params.outSAMtype_sort = 'Unsorted'
star_align_args       += ' --outSAMtype ' + params.outSAMtype_file + ' ' + params.outSAMtype_sort + ' '


/* ========================================================================================
    SUBREAD (FEATURECOUNTS) PARAMETERS
======================================================================================== */
// B flag
params.featurecounts_B_flag = true 
// Only count read pairs that have both ends aligned.

if(params.featurecounts_B_flag){
    featurecounts_args += " -B "
}

// C flag
params.featurecounts_C_flag = true
// -C  Do not count read pairs that have their two ends mapping
//     to different chromosomes or mapping to same chromosome
//     but on different strands.

if(params.featurecounts_C_flag){
    featurecounts_args += " -C "
}

// strandness
if (params.strandness == 'forward') {
    featurecounts_args += " -s 1 "
} else if (params.strandness == 'reverse') {
    featurecounts_args += " -s 2 "
} else if (params.strandness == 'unstranded' || params.strandness == 'smartseq2') {
    featurecounts_args += " -s 0 "
}


/* ========================================================================================
    GENOMES
======================================================================================== */
params.custom_genome_file = '' // Option to add a directory for a custom genome file

include { getGenome } from './modules/genomes.mod.nf' params(custom_genome_file: params.custom_genome_file)
genome = getGenome(params.genome)


/* ========================================================================================
    FILES CHANNEL
======================================================================================== */
include { makeFilesChannel; getFileBaseNames } from './modules/files.mod.nf'
file_ch = makeFilesChannel(args)


/* ========================================================================================
    MESSAGES
======================================================================================== */
// Validate strandedness
assert params.strandness == 'forward' || params.strandness == 'reverse' || params.strandness == 'unstranded' || params.strandness == 'smartseq2' : "Invalid strand orientation option: >>${params.strandness}<<. Valid options are: 'forward', 'reverse', 'unstranded' or 'smartseq2'\n\n"
println ("Using strand orientation: " + params.strandness)

// Validate aligner
assert params.aligner == 'star' || params.aligner == 'hisat2' : "Invalid aligner option: >>${params.aligner}<<. Valid options are: 'star' or 'hisat2'\n\n"
println ("Using aligner: " + params.aligner)


/* ========================================================================================
    WORKFLOW
======================================================================================== */
include { FASTQC }                     from './modules/fastqc.mod.nf'
include { FASTQC as FASTQC2 }          from './modules/fastqc.mod.nf'
include { FASTQ_SCREEN }               from './modules/fastq_screen.mod.nf' params(fastq_screen_conf: params.fastq_screen_conf)
include { TRIM_GALORE }                from './modules/trim_galore.mod.nf'
include { HISAT2 }                     from './modules/hisat2.mod.nf'       params(genome: genome, bam_output: false)
include { STAR_ALIGN }                 from './modules/star.mod.nf'         params(genome: genome, bam_output: false)
include { SAMTOOLS_VIEW }              from './modules/samtools.mod.nf'
include { SAMTOOLS_SORT }              from './modules/samtools.mod.nf'
include { SAMTOOLS_INDEX }             from './modules/samtools.mod.nf'
include { FEATURECOUNTS }              from './modules/subread.mod.nf'      params(genome: genome)
include { FEATURECOUNTS_MERGE_COUNTS } from './modules/subread.mod.nf'
include { MULTIQC }                    from './modules/multiqc.mod.nf' 

workflow {

    main:

        // STAR_ALIGN aligner
        if (params.aligner == 'star'){

            if (!params.skip_qc){ 

                FASTQC                          (file_ch, outdir, params.fastqc_args)
                TRIM_GALORE                     (file_ch, outdir, params.trim_galore_args)

                if (!params.skip_fastq_screen){ 
                FASTQ_SCREEN                    (TRIM_GALORE.out.reads, outdir, params.fastq_screen_args)
                }

                FASTQC2                         (TRIM_GALORE.out.reads, outdir, params.fastqc_args)
                STAR_ALIGN                      (TRIM_GALORE.out.reads, outdir, params.star_align_args)

            } else {
                STAR_ALIGN                      (TRIM_GALORE.out.reads, outdir, params.star_align_args)
            }
        
            SAMTOOLS_VIEW               (STAR_ALIGN.out.bam, outdir, params.samtools_view_args)
            SAMTOOLS_SORT               (SAMTOOLS_VIEW.out.bam, outdir, params.samtools_sort_args)
            SAMTOOLS_INDEX              (SAMTOOLS_SORT.out.bam, outdir, params.samtools_index_args)
            FEATURECOUNTS               (SAMTOOLS_SORT.out.bam, STAR_ALIGN.out.single_end, outdir, params.featurecounts_args)
            featurecounts_merge_counts_ch = FEATURECOUNTS.out.counts.collect()
            FEATURECOUNTS_MERGE_COUNTS  (featurecounts_merge_counts_ch, outdir)
        }

        // HISAT2 aligner
        if (params.aligner == 'hisat2'){

            if (!params.skip_qc){ 

                FASTQC                          (file_ch, outdir, params.fastqc_args)
                TRIM_GALORE                     (file_ch, outdir, params.trim_galore_args)

                if (!params.skip_fastq_screen){ 
                FASTQ_SCREEN                    (TRIM_GALORE.out.reads, outdir, params.fastq_screen_args)
                }

                FASTQC2                         (TRIM_GALORE.out.reads, outdir, params.fastqc_args)
                HISAT2                          (TRIM_GALORE.out.reads, outdir, params.hisat2_args)

            } else {
                HISAT2                          (TRIM_GALORE.out.reads, outdir, params.hisat2_args)
            }
        
            SAMTOOLS_VIEW               (HISAT2.out.bam, outdir, params.samtools_view_args)
            SAMTOOLS_SORT               (SAMTOOLS_VIEW.out.bam, outdir, params.samtools_sort_args)
            SAMTOOLS_INDEX              (SAMTOOLS_SORT.out.bam, outdir, params.samtools_index_args)
            FEATURECOUNTS               (SAMTOOLS_SORT.out.bam, HISAT2.out.single_end, outdir, params.featurecounts_args)
            featurecounts_merge_counts_ch = FEATURECOUNTS.out.counts.collect()
            FEATURECOUNTS_MERGE_COUNTS  (featurecounts_merge_counts_ch, outdir)
        }


        /* ========================================================================================
            Reports
        ======================================================================================== */

        // Merging channels for MultiQC
        if (!params.skip_qc){

            multiqc_ch = FASTQC.out.report.mix(
                         TRIM_GALORE.out.report,
                         FASTQC2.out.report.ifEmpty([])
                         ).collect()

            if (!params.skip_fastq_screen){
                multiqc_ch = multiqc_ch.mix(
                            FASTQ_SCREEN.out.report.ifEmpty([])
                            ).collect()
            }

        } else {
            if (params.aligner == 'star'){
                multiqc_ch = STAR_ALIGN.out.stats.ifEmpty([])
            }
            if (params.aligner == 'hisat2'){
                multiqc_ch = HISAT2.out.stats.ifEmpty([])
            }
        }

        multiqc_ch = multiqc_ch.mix(
                    FEATURECOUNTS.out.summary.ifEmpty([])
                    ).collect() 

        MULTIQC (multiqc_ch, outdir, params.multiqc_args)
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Jobname     : ${workflow.runName}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
    .stripIndent()

    sendMail(to: "${workflow.userName}@ethz.ch", subject: 'Minimal pipeline execution report', body: msg)
}
