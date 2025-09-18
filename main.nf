#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    INPUT FILES
======================================================================================== */
params.input = null

if (!params.input) {
    error "Input not specified. Use --input to specify the input."
}

input_files = file(params.input)


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
params.skip_qc             = false
params.skip_fastq_screen   = false
params.skip_quantification = false


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.fastq_screen_conf = '/cluster/work/nme/software/config/fastq_screen.conf' // FastQ Screen config file directory
params.genome            = 'Mus_musculus_GRCm39' // Default genome
params.strandness        = 'smartseq2' // Library strandedness
params.single_end        = false // Force to input files to be single-end


/* ========================================================================================
    EXTRA PARAMETERS
======================================================================================== */
params.fastqc_args         = ''
params.fastq_screen_args   = ''
params.trim_galore_args    = ''
params.star_align_args     = ''
params.hisat2_align_args   = ''
params.salmon_quant_args   = ''
params.featurecounts_args  = ''
params.multiqc_args        = ''
params.samtools_sort_args  = ''
params.samtools_index_args = ''

fastqc_args         = params.fastqc_args       
fastq_screen_args   = params.fastq_screen_args
trim_galore_args    = params.trim_galore_args 
star_align_args     = params.star_align_args 
hisat2_align_args   = params.hisat2_align_args 
salmon_quant_args   = params.salmon_quant_args
featurecounts_args  = params.featurecounts_args
multiqc_args        = params.multiqc_args
samtools_sort_args  = params.samtools_sort_args
samtools_index_args = params.samtools_index_args


/* ========================================================================================
    ALIGNER
======================================================================================== */
params.aligner = 'star'


/* ========================================================================================
    HISAT2_ALIGN PARAMETERS
======================================================================================== */
// no-softclip: no soft-clipping
params.hisat2_no_softclip = true

if(params.hisat2_no_softclip){
    hisat2_align_args += " --no-softclip "
}

// no-mixed: suppress unpaired alignments for paired reads
params.hisat2_no_mixed = true

if(params.hisat2_no_mixed){
    hisat2_align_args += " --no-mixed "
}

// no-discordant: suppress discordant alignments for paired reads
params.hisat2_no_discordant = true

if(params.hisat2_no_discordant){
    hisat2_align_args += " --no-discordant "
}


/* ========================================================================================
    STAR PARAMETERS
======================================================================================== */
// outSAMtype
params.star_outSAMtype_file = 'BAM'
params.star_outSAMtype_sort = 'Unsorted'
star_align_args       += ' --outSAMtype ' + params.star_outSAMtype_file + ' ' + params.star_outSAMtype_sort + ' '


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
file_ch = makeFilesChannel(input_files)


/* ========================================================================================
    MESSAGES
======================================================================================== */
// Validate strandedness
assert params.strandness == 'forward' || params.strandness == 'reverse' || params.strandness == 'unstranded' || params.strandness == 'smartseq2' : "Invalid strand orientation option: >>${params.strandness}<<. Valid options are: 'forward', 'reverse', 'unstranded' or 'smartseq2'\n\n"
println ("Using strand orientation: " + params.strandness)

// Validate aligner
assert params.aligner == 'star' || params.aligner == 'hisat2' || params.aligner == 'salmon' : "Invalid aligner option: >>${params.aligner}<<. Valid options are: 'star' or 'hisat2' or 'salmon'\n\n"
println ("Using aligner: " + params.aligner)


/* ========================================================================================
    WORKFLOW
======================================================================================== */
include { FASTQC }                               from './modules/fastqc.mod.nf'
include { FASTQC as FASTQC2 }                    from './modules/fastqc.mod.nf'
include { FASTQ_SCREEN }                         from './modules/fastq_screen.mod.nf' params(fastq_screen_conf: params.fastq_screen_conf)
include { TRIM_GALORE }                          from './modules/trim_galore.mod.nf'
include { HISAT2_ALIGN }                         from './modules/hisat2.mod.nf'       params(genome: genome, bam_output: false)
include { STAR_ALIGN }                           from './modules/star.mod.nf'         params(genome: genome, bam_output: false)
include { SALMON_QUANT }                         from './modules/salmon.mod.nf'       params(genome: genome, bam_output: false)
include { SAMTOOLS_SORT }                        from './modules/samtools.mod.nf'
include { SAMTOOLS_INDEX }                       from './modules/samtools.mod.nf'
include { FEATURECOUNTS }                        from './modules/subread.mod.nf'      params(genome: genome)
include { FEATURECOUNTS_MERGE_COUNTS }           from './modules/subread.mod.nf'
include { FEATURECOUNTS_MERGE_COUNTS_SALMON }    from './modules/subread.mod.nf'
include { FEATURECOUNTS_MERGE_COUNTS_SALMON_TX } from './modules/subread.mod.nf'
include { MULTIQC }                              from './modules/multiqc.mod.nf' 

workflow {

    main:

        // STAR_ALIGN aligner
        if (params.aligner == 'star'){

            if (!params.skip_qc){ 

                FASTQC                          (file_ch, outdir, fastqc_args)

                if (!params.skip_fastq_screen){ 
                FASTQ_SCREEN                    (file_ch, outdir, fastq_screen_args)
                }

                TRIM_GALORE                     (file_ch, outdir, trim_galore_args)
                FASTQC2                         (TRIM_GALORE.out.reads, outdir, fastqc_args)
                STAR_ALIGN                      (TRIM_GALORE.out.reads, outdir, star_align_args)

            } else {
                STAR_ALIGN                      (file_ch, outdir, star_align_args)
            }
        
            SAMTOOLS_SORT               (STAR_ALIGN.out.bam, outdir, samtools_sort_args)
            SAMTOOLS_INDEX              (SAMTOOLS_SORT.out.bam, outdir, samtools_index_args)

            if (!params.skip_quantification){
                FEATURECOUNTS               (SAMTOOLS_SORT.out.bam, STAR_ALIGN.out.single_end, outdir, featurecounts_args)
                featurecounts_merge_counts_ch = FEATURECOUNTS.out.counts.collect()
                FEATURECOUNTS_MERGE_COUNTS  (featurecounts_merge_counts_ch, outdir)
            }
        }

        // HISAT2_ALIGN aligner
        if (params.aligner == 'hisat2'){

            if (!params.skip_qc){ 

                FASTQC                          (file_ch, outdir, fastqc_args)
                
                if (!params.skip_fastq_screen){ 
                FASTQ_SCREEN                    (file_ch, outdir, fastq_screen_args)
                }

                TRIM_GALORE                     (file_ch, outdir, trim_galore_args)
                FASTQC2                         (TRIM_GALORE.out.reads, outdir, fastqc_args)
                HISAT2_ALIGN                    (TRIM_GALORE.out.reads, outdir, hisat2_align_args)

            } else {
                HISAT2_ALIGN                    (file_ch, outdir, hisat2_align_args)
            }

            SAMTOOLS_SORT               (HISAT2_ALIGN.out.bam, outdir, samtools_sort_args)
            SAMTOOLS_INDEX              (SAMTOOLS_SORT.out.bam, outdir, samtools_index_args)

            if (!params.skip_quantification){
                FEATURECOUNTS               (SAMTOOLS_SORT.out.bam, HISAT2_ALIGN.out.single_end, outdir, featurecounts_args)
                featurecounts_merge_counts_ch = FEATURECOUNTS.out.counts.collect()
                FEATURECOUNTS_MERGE_COUNTS  (featurecounts_merge_counts_ch, outdir)
            }
        }


        // SALMON_QUANT aligner
        if (params.aligner == 'salmon'){

            if (!params.skip_qc){ 

                FASTQC                          (file_ch, outdir, fastqc_args)
                
                if (!params.skip_fastq_screen){ 
                FASTQ_SCREEN                    (file_ch, outdir, fastq_screen_args)
                }

                TRIM_GALORE                     (file_ch, outdir, trim_galore_args)
                FASTQC2                         (TRIM_GALORE.out.reads, outdir, fastqc_args)
                SALMON_QUANT                    (TRIM_GALORE.out.reads, outdir, salmon_quant_args, params.strandness)

            } else {
                SALMON_QUANT                    (file_ch, outdir, salmon_quant_args, params.strandness)
            }
            if (!params.skip_quantification){
                featurecounts_merge_counts_ch         = SALMON_QUANT.out.counts_gene.collect()
                FEATURECOUNTS_MERGE_COUNTS_SALMON     (featurecounts_merge_counts_ch, outdir)
                featurecounts_merge_counts_tx_ch      = SALMON_QUANT.out.counts_tx.collect()
                FEATURECOUNTS_MERGE_COUNTS_SALMON_TX  (featurecounts_merge_counts_tx_ch, outdir)

            }

        }


        /* ========================================================================================
            Reports
        ======================================================================================== */

        // Merging channels for MultiQC
        if (!params.skip_qc){

            multiqc_ch = FASTQC.out.report.mix(
                        TRIM_GALORE.out.report.ifEmpty([]),
                        FASTQC2.out.report.ifEmpty([])
                        ).collect()

            if (!params.skip_fastq_screen){
                multiqc_ch = multiqc_ch.mix(
                            FASTQ_SCREEN.out.report.ifEmpty([])
                            ).collect()
            }

            if (params.aligner == 'star'){
                multiqc_ch = multiqc_ch.mix(
                            STAR_ALIGN.out.log_final.ifEmpty([])
                            ).collect() 
            }
            if (params.aligner == 'hisat2'){
                multiqc_ch = multiqc_ch.mix(
                            HISAT2_ALIGN.out.stats.ifEmpty([])
                            ).collect()
            }

        } else {

            if (params.aligner == 'star'){
                multiqc_ch = STAR_ALIGN.out.log_final.ifEmpty([]).collect() 
            }
            if (params.aligner == 'hisat2'){
                multiqc_ch = HISAT2_ALIGN.out.stats.ifEmpty([]).collect()
            }

        }

        if (!params.skip_quantification){
            multiqc_ch = multiqc_ch.mix(
                        FEATURECOUNTS.out.summary.ifEmpty([])
                        ).collect()
        }

        //MULTIQC (multiqc_ch, outdir, multiqc_args)
}
