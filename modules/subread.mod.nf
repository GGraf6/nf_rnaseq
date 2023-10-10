#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* ========================================================================================
    PROCESSES
======================================================================================== */

// FEATURECOUNTS
process FEATURECOUNTS {	

    label "featureCounts"
	tag "$bam" // Adds name to job submission

	input:
		path(bam)
        val(single_end)
		val(outputdir)
		val(featurecounts_args)

	output:
        path("*featureCounts.txt")        , emit: counts
        path("*featureCounts.txt.summary"), emit: summary
        
	 	publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true

	script:
		/* ==========
			Paired-end or single-end
		========== */
        paired_end = single_end ? '' : '-p'

		/* ==========
			Basename
		========== */
        basename = bam.toString() - ".bam"

		/* ==========
			Annotation file
		========== */
        annotation = params.genome["gtf"]
        
		"""
		module load subread

		featureCounts \\
            $featurecounts_args \\
            $paired_end \\
            -T $task.cpus \\
            -a $annotation \\
            -o ${basename}.featureCounts.txt \\
            ${bam}
		"""
}

// FEATURECOUNTS_MERGE_COUNTS
process FEATURECOUNTS_MERGE_COUNTS {

    input:
        path('counts/*')
        val(outputdir)

    output:
        path "*.txt", emit: merged_counts

        publishDir "$outputdir/aligned/counts", mode: "link", overwrite: true

    script:
    """
    mkdir -p tmp/counts

    cut -f 1 `ls ./counts/* | head -n 1` | grep -v "^#" > ids.tsv

    for fileid in `ls ./counts/*featureCounts.txt`; do
        samplename=`basename \$fileid | sed s/\\.featureCounts.txt\$//g`
        echo \$samplename > tmp/counts/\$samplename.featureCounts.txt
        grep -v "^#" \${fileid} | cut -f 7 | tail -n+2 >> tmp/counts/\$samplename.featureCounts.txt
    done
    
    paste ids.tsv tmp/counts/* > gene_counts_merged.txt
    """
}
