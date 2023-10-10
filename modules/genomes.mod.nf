#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PARAMETERS
======================================================================================== */
params.custom_genome_file = ''


/* ========================================================================================
    FUNCTIONS
======================================================================================== */

// getGenome
def getGenome(name) {

    // Find a file with the same name as the genome in our genomes directory
    scriptDir = workflow.projectDir

    if (params.custom_genome_file){
        
        fileName = params.custom_genome_file

    } else {

        fileName = scriptDir.toString() + "/genomes/" + name + ".genome"

    }

    // die gracefully if the user specified an incorrect genome
    def testFile = new File(fileName)
    if (!testFile.exists()) {
        println("\nFile >>$fileName<< does not exist. Listing available genomes...\n")
        listGenomes()
    }
    else {
        // println ("File $fileName exists.")
    }

    genomeFH = new File (fileName).newInputStream()

    genomeValues = [:]  // initialising map. name is also part of each .genome file

    genomeFH.eachLine {
        sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}


// listGenomes
def listGenomes() {

    println ("These genomes are currently available to choose from:")
    println ("=====================================================")
    scriptDir = workflow.projectDir + "/genomes/"

    allFiles = scriptDir.list()

    for( def file : allFiles.sort() ) {

        if( file =~ /.genome$/){

            genomeFH = new File(scriptDir.toString() + "/$file").newInputStream()
            name = file.replaceFirst(/.genome/, "")

            println (name)
            genomeFH.eachLine {
                if (params.verbose){
                    println ("\t$it")
                }
            }
        }
    }
    println ("\nTo see this list of available genomes with more detailed information about paths and indexes,\nplease re-run the command including '--list_genomes --verbose'\n\n")

    System.exit(1)
}
