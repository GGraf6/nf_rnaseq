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
			File names
		========== */
		readString = ""
		if (reads instanceof List) {
			readString  = "-1 " + reads[0] + " -2 " + reads[1]
			single_end  = false
		}
		else {
			readString  = "-r " + reads
			single_end  = true
		}

        /* ==========
			Index
		========== */
		index = params.genome["salmon"]


		/* ==========
			Transcript to gene ID conversion
		========== */
		anno = params.genome["tx_to_gn"]


// TODO: 1) provide gtf or tx_to_gn? 2) parse strand based on PE/SE and read orientation 3) output?

		"""
		salmon quant -l $strand --threads ${task.cpus} --geneMap $gtf -i $index ${readString} -o $outfile
        //salmon quant --geneMap ${anno} --threads ${task.cpus} --libType=$strandedness $reference ${reads} ${salmon_quant_args} -o ${basename}.salmon.txt
		"""
}
