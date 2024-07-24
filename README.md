# RNA Sequencing Pipeline

<img width="30%" src="https://raw.githubusercontent.com/nextflow-io/trademark/master/nextflow-logo-bg-light.png" />
<img width="30%" src="https://tower.nf/assets/nf-tower-black.svg" />

A Nextflow pipeline to perform quality control, alignment, and quantification of RNA sequencing data.

>The pipeline was created to run on the [ETH Euler cluster](https://scicomp.ethz.ch/wiki/Euler) and it relies on the server's genome files. Thus, the pipeline needs to be adapted before running it in a different HPC cluster.

## Pipeline steps
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
3. [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
4. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
5. [STAR](https://github.com/alexdobin/STAR) or [HISAT2](https://daehwankimlab.github.io/hisat2/)
6. [Samtools sort](https://www.htslib.org/doc/samtools-sort.html)
7. [Samtools index](https://www.htslib.org/doc/samtools-index.html)
8. [featureCounts](https://subread.sourceforge.net/featureCounts.html#:~:text=featureCounts%20is%20a%20highly%20efficient,and%20genomic%20DNA%2Dseq%20reads.)
19. [MultiQC](https://multiqc.info/)

## Required parameters

Path to the folder where the FASTQ files are located.
``` bash
--input /cluster/work/nme/data/josousa/project/fastq/*fastq.gz
```

Output directory where the files will be saved.
``` bash
--outdir /cluster/work/nme/data/josousa/project
```

### Input optional parameters

- Option to force the pipeline to assign input as single-end.

    `--single_end`

    >_By default, the pipeline detects whether the input files are single-end or paired-end._


- Option to select RNA-Seq library strandness. This will only affect quantification.
    ``` bash
    --strandness 'smartseq2' # Default (same as 'unstranded')
    --strandness 'forward'
    --strandness 'reverse'
    --strandness 'unstranded'
    ```

    _This option will only affect quantification._

## Genomes
- Reference genome used for alignment.

    `--genome`

    Available genomes:
    ``` bash
        GRCm39 # Default
        GRCm38
        GRCh38
        GRCh37 
        panTro6
        CHIMP2.1.4
        BDGP6
        susScr11
        Rnor_6.0
        R64-1-1
        TAIR10
        WBcel235
        E_coli_K_12_DH10B
        E_coli_K_12_MG1655
        Vectors
        Lambda
        PhiX
        Mitochondria
    ```

- Option to use a custom genome for alignment by providing an absolute path to a custom genome file.

    ``` bash
    --custom_genome_file '/cluster/work/nme/data/josousa/project/genome/CHM13.genome'
    ```

    Example of a genome file:
    ``` bash    
    name           GRCm39
    species        Mouse
    star           /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/STARIndex/
    hisat2         /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/genome
    hisat2_splices /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Sequence/Hisat2Index/splice_sites.txt
    gtf            /cluster/work/nme/genomes/Mus_musculus/Ensembl/GRCm39/Annotation/Genes/genes.gtf
    ```

## Aligner options

- Option to choose the aligner.
    ``` bash
    --aligner 'star' # Default
    --aligner 'hisat2'
    ```

### HISAT2 parameters
- Option to choose no soft-clipping.

    `--hisat2_no_softclip` _Default: true_

- Option to suppress unpaired alignments for paired reads

    `--hisat2_no_mixed` _Default: true_

- Option to suppress discordant alignments for paired reads.

    `--hisat2_no_discordant` _Default: true_


## FastQ Screen optional parameters

- Option to provide a custom FastQ Screen config file.
    ``` bash
    --fastq_screen_conf '/cluster/work/nme/software/config/fastq_screen.conf' # Default
    ```

- Option to pass the flag --bisulfite to FastQ Screen.

    `--bisulfite` _Default: false_

## featureCounts optional parameters

- Option to only count read pairs that have both ends aligned.

    `--featurecounts_B_flag` _Default: true_

- Option to not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.

    `--featurecounts_C_flag` _Default: true_


## Skipping options
- Option to skip FastQC, TrimGalore, and FastQ Screen. The first step of the pipeline will be the Bismark alignment. 
`--skip_qc`

- Option to skip FastQ Screen. 
`--skip_fastq_screen`

- Option to skip quantification. 
`--skip_quantification`

## Extra arguments

- Option to add extra arguments to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
`--fastqc_args`

- Option to add extra arguments to [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/).
`--fastq_screen_args`

- Option to add extra arguments to [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
`--trim_galore_args`

- Option to add extra arguments to the [STAR](https://github.com/alexdobin/STAR) aligner.
`--star_align_args`

- Option to add extra arguments to the [HISAT2](https://daehwankimlab.github.io/hisat2/) aligner.
`--hisat2_align_args`

- Option to add extra arguments to [Samtools sort](https://www.htslib.org/doc/samtools-sort.html).
`--samtools_sort_args`

- Option to add extra arguments to [Samtools index](https://www.htslib.org/doc/samtools-index.html).
`--samtools_index_args`

- Option to add extra arguments to [featureCounts](https://subread.sourceforge.net/featureCounts.html#:~:text=featureCounts%20is%20a%20highly%20efficient,and%20genomic%20DNA%2Dseq%20reads.).
`--featurecounts_args`

- Option to add extra arguments to [MultiQC](https://multiqc.info/).
`--multiqc_args`

## Acknowledgements
This pipeline was adapted from the Nextflow pipelines created by the [Babraham Institute Bioinformatics Group](https://github.com/s-andrews/nextflow_pipelines) and from the [nf-core](https://nf-co.re/) pipelines. We thank all the contributors for both projects. We also thank the [Nextflow community](https://nextflow.slack.com/join) and the [nf-core community](https://nf-co.re/join) for all the help and support.

