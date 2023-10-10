#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.bam_output = true // Setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */

// SAMTOOLS_VIEW
process SAMTOOLS_VIEW {	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission

	input:
		path(bam)
		val(outputdir)
		val(samtools_view_args)

	output:
		path "*bam", emit: bam

    script:
		/* ==========
			Basename
		========== */
		basename = bam.toString() - ".bam"

		"""
		module load samtools

		samtools view --threads ${task.cpus} $samtools_view_args -bS -F 4 -F 8 -F 256 $bam -o ${basename}.bam
		"""
}

// SAMTOOLS_SORT
process SAMTOOLS_SORT {	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission

	input:
		path(bam)
		val(outputdir)
		val(samtools_sort_args)

	output:
		path "*bam", emit: bam
		publishDir "$outputdir/aligned/bam", mode: "link", overwrite: true, enabled: params.bam_output

    script:
		/* ==========
			Basename
		========== */
		basename = bam.toString() - ".bam"

		"""
		module load samtools

		samtools sort --threads ${task.cpus} $samtools_sort_args $bam -o ${basename}.sorted.bam
    	"""
}

// SAMTOOLS_INDEX
process SAMTOOLS_INDEX {	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission

	input:
		path(bam)
		val(outputdir)
		val(samtools_index_args)

	output:
		path "*.bai", emit: bai
		publishDir "$outputdir/aligned/bam", mode: "link", overwrite: true

    script:

		"""
		module load samtools

		samtools index -@ ${task.cpus} $samtools_index_args $bam
		"""
}
