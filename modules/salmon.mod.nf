#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */

process SALMON_QUANT {

	label 'salmon'
	tag "$name" // Adds name to job submission

	container 'docker://combinelab/salmon:1.10.3'

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(salmon_quant_args)

	output:
		path("${prefix}") , emit: results
		path("*info.json"), emit: json_info, optional: true

		publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true

	script:

		/* ==========
			GTF file
		========== */
		gtf = params.genome["gtf"]

		"""
		salmon quant --geneMap ${gtf} --threads ${task.cpus} --libType=$strandedness $reference ${reads} ${salmon_quant_args} -o ${basename}.salmon.txt
		"""
}
