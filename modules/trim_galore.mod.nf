#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.clip_r2 = false


/* ========================================================================================
    PROCESSES
======================================================================================== */
process TRIM_GALORE {

	label 'trimGalore'
	tag "$name" // Adds name to job submission

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(trim_galore_args)

	output:
		tuple val(name), path("*fq.gz"), 			 emit: reads
		path "*trimming_report.txt", optional: true, emit: report
			
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*fq.gz"
		publishDir "$outputdir/unaligned/logs",  mode: "link", overwrite: true, pattern: "*trimming_report.txt"

	script:
	
		/* ==========
			Paired-end
		========== */
		if (reads instanceof List){
			trim_galore_args += " --paired "

			if (params.clip_r2){
				trim_galore_args += " --clip_r2 ${params.clip_r2} "
			}
		}

		"""
		module load trimgalore

		trim_galore --cores 4 $trim_galore_args ${reads}
		"""
}
