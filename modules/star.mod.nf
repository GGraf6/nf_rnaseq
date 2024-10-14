#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published

params.gzip       = false // To detect if the input files are gzip compressed
gzip              = params.gzip


/* ========================================================================================
    PROCESSES
======================================================================================== */
process STAR_ALIGN {

	label 'star_align'
	tag "$name" // Adds name to job submission

	container 'docker://josousa/star:2.5.1b'

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(star_align_args)

	output:
		path('*d.out.bam')       , emit: bam
		path('*Log.final.out')   , emit: log_final
		path('*Log.out')         , emit: log_out
		path('*Log.progress.out'), emit: log_progress

		path('*sortedByCoord.out.bam')  , optional: true, emit: bam_sorted
		path('*toTranscriptome.out.bam'), optional: true, emit: bam_transcript
		path('*Aligned.unsort.out.bam') , optional: true, emit: bam_unsorted
		path('*fastq.gz')               , optional: true, emit: fastq
		path('*.tab')                   , optional: true, emit: tab
		path('*.out.junction')          , optional: true, emit: junction
		path('*.out.sam')               , optional: true, emit: sam

		val(single_end), emit: single_end
		
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*out"

	script:
		/* ==========
			File names
		========== */
		readString = ""
		if (reads instanceof List) {
			readString  = reads[0] + " " + reads[1]
			single_end  = false

			if (reads[0].getExtension() == 'gz'){
				gzip = true
			}
		}
		else {
			readString  = reads
			single_end  = true

			if (reads.getExtension() == 'gz'){
				gzip = true
			}
		}

		/* ==========
			gzip input files
		========== */
		if (gzip){
			star_align_args += " --readFilesCommand zcat "
		}

		/* ==========
			Index
		========== */
		index = params.genome["star"]
		
		/* ==========
			Basename
		========== */
		star_name = name + "_" + params.genome["name"] + "_" + "star" +  "."

		"""
		STAR \\
			--genomeDir ${index} \\
			--outFileNamePrefix ${star_name} \\
			--runThreadN ${task.cpus} \\
			${star_align_args} \\
			--readFilesIn ${readString}
		"""
}
