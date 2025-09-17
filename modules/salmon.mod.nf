#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */

process SALMON_QUANT {

	label 'salmon_quant'
	tag "$name" // Adds name to job submission

	//container 'docker://combinelab/salmon:1.10.3'

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(salmon_quant_args)
        val(strandness)

	output:
		path('*/quant.genes.sf')    , emit: counts_gene
		path('*/quant.sf')          , emit: counts_tx
		path('*/*json')             , emit: jsons

		publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true, pattern: "*/*sf"
        publishDir "$outputdir/aligned/logs",   mode: "link", overwrite: true, pattern: "*/*json"

	script:

        /* ==========
			File names and strandness
		========== */
		readString   = ""
        strandString = ""
		if (reads instanceof List) {
			readString   = "-1 " + reads[0] + " -2 " + reads[1]
            strandString += "I"
		}
		else {
			readString   = "-r " + reads
		}


        /* ==========
			strandness: this depends also on the read SE or PE configuration
		========== */
        if (strandness == 'forward') {
            strandString += "SF"
        } else if (strandness == 'reverse') {
            strandString += "SR"
        } else if (strandness == 'unstranded' || params.strandness == 'smartseq2') {
            strandString += "U"
        }


        /* ==========
			Index
		========== */
        index = params.genome["salmon"]


		/* ==========
			Transcript to gene ID conversion
		========== */
		tx_to_gene = params.genome["tx_to_gene"]


        /* ==========
			Basename
		========== */
		salmon_name = name + "_" + params.genome["name"] + "_" + "salmon"



		"""
		module load salmon
		salmon quant -l ${strandString} --threads ${task.cpus} --geneMap ${tx_to_gene} -i ${index} ${readString} ${salmon_quant_args} -o ./${salmon_name}
		# rename salmon outputs
		mv \${salmon_name}/quant.sf \${salmon_name}/\${salmon_name}_quant.sf
		mv \${salmon_name}/quant.genes.sf \${salmon_name}/\${salmon_name}_quant.genes.sf
		"""
}
