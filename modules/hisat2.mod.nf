#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */
process HISAT2 {

	label 'hisat2'
	tag "$name" // Adds name to job submission

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(hisat2_args)

	output:
		path "*bam",       emit: bam
		path "*stats.txt", emit: stats
		val(single_end)  , emit: single_end
		
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*stats.txt"

	script:
		/* ==========
			File names
		========== */
		readString = ""
		if (reads instanceof List) {
			readString  = "-1 " + reads[0] + " -2 " + reads[1]
			hisat2_args = hisat2_args + " --no-mixed --no-discordant"
			single_end  = false
		}
		else {
			readString  = "-U " + reads
			single_end  = true
		}

		/* ==========
			Index
		========== */
		index = params.genome["hisat2"]
		
		/* ==========
			Splices
		========== */
		// TODO: need to add a check if the splice-site infile is present or not, and leave out this parameter otherwise
		splices = " --known-splicesite-infile " + params.genome["hisat2_splices"]

		/* ==========
			Basename
		========== */
		hisat_name = name + "_" + params.genome["name"]

		"""
		module load hisat2 samtools

		hisat2 -p ${task.cpus} ${hisat2_args} -x ${index} ${splices} ${readString} 2>${hisat_name}_ht2_stats.txt -> ${hisat_name}_ht2.bam
		"""
}
