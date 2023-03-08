#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
* Pipeline info and input parameter paths
*/

log.info """\

    P R O T O N F _ n o v o PIPELINE
    ========================
    followup to protoNF1 prototype that incorporates the novoalign aligner package
    ========================
    Reads                   : ${params.reads}
    Barcodes                : ${params.barcodes}
    Novoalign Index         : ${params.novoindex}
    GTF Annotation file     : ${params.gtf}
    Genome Fasta            : ${params.genome}

    """
    .stripIndent()

/*
* Adapter Trimming - Flexbar to remove the Illumina or other sequencing adapters from the end of the reads
*/

process FLEXBAR {
    input: 
    path(reads)

    output: 
    file('flexbar_trimmed.fastq')

    script:
    """
    flexbar -r $reads -qf i1.8 -n 10 -ao 7 --output-reads flexbar_trimmed.fastq -qt 30
    """
} 

/*
* Demultiplex - pyBarcodeFilter to split to the sequence data into the respective experiments they were run in
*/

process DEMULTIPLEX {

    input:
    path(flexbar_output)
    path(barcodes)

    output:
    path "flexbar_trimmed_*_*.fastq"

    script:
    """
    pyBarcodeFilter.py -f '$flexbar_output' -b '$barcodes' -m '1'
    """
    //run pyBarcode Filter on the trimmed reads from previous step using the barcodes provides with a maximum mismatch of 1 (this is the default)
} 

/*
* Collapse - pyFastDuplicateRemover to remove duplicated while retaining the sequence variation
*/

process COLLAPSE{

    input:
    path(fastq_files)

    output:
    path "flexbar_trimmed_*_collapsed.fasta"

    script:
    """
    pyFastqDuplicateRemover.py -f '$fastq_files' -o '${fastq_files.baseName}_collapsed.fasta'
    """
} 

/*
*  Test Novoalign, aligns preprocessed files to novoindex, also pregenerated and fed as an input
*/

process NOVOALIGN {

    publishDir './results/Alignment' , mode: 'copy', overwrite: false

    input:
    path(novoindex)
    each(collapsed_files)

    output:
    path "*.novo"

    script:
    """
    novoalign -d '$novoindex' -f '$collapsed_files' -r Random > '${collapsed_files.baseName}_aligned.novo'
    """
}

/*
* pyReadCounter, generate hit tables, show number of reads mapped to reference and what theyre mapped to 
*/

process RUNPYREADCOUNTERS{

    publishDir './results/hitTable' , mode: 'copy', overwrite: false

    input:
    path(gtf)
    each(novoaligned_files)

    output:
    path "*.gtf"

    script:
    """
    pyReadCounters.py -f '$novoaligned_files' --gtf '$gtf' -v --rpkm -o '${novoaligned_files.baseName}' 
    """
} 

process RUNSECONDPYREADCOUNTERS{

    publishDir './results/hitTable2' , mode: 'copy', overwrite: false

    input:
    path(gtf)
    each(novoaligned_files)

    output:
    path "*.gtf"

    script:
    """
    pyReadCounters.py -f '$novoaligned_files' --gtf '$gtf' -v --rpkm -o '${novoaligned_files.baseName}_nomuts' --mutations  nomuts --blocks
    """
} 

/*
* pyCalculateChromosomeLengths.py, generation of a chromosome file for usage by the subsequent 2 steps
*/

process CHROMOSOMELENGTH{

    input:
    path(genome)

    output:
    path "*.txt"

    script:
    """
    pyCalculateChromosomeLengths.py -f '$genome' -o '${genome.baseName}_chromosome.txt' --file_type=fasta
    """
}

/*
* pyGTF2sgr.py - takes from the pyReadCounters GTF file and generates coverage sgr file 
*/

process COVERAGE{

    publishDir './results/sgrCoverage' , mode: 'copy', overwrite: false

    input:
    path(chromosome_files)
    each(secondPyReadCounter_files)
    
    

    output:
    file "*.sgr" 

    script:
    """
    pyGTF2sgr.py --gtf '$secondPyReadCounter_files' --zeros --count -v -o '${secondPyReadCounter_files.baseName}' -c '$chromosome_files'
    """

}

/*
* takes pyReadCounters GTF file and generates a readable Bedgraph file
*/

process BEDGRAPH{

    publishDir './results/bedgraph' , mode: 'copy', overwrite: false

    input:
    path(chromosome_files)
    each(secondPyReadCounter_files)
    
    

    output:
    path "*.bedgraph"

    script:
    """
    pyGTF2bedGraph.py --gtf '$secondPyReadCounter_files' --count -v --permillion -o '${secondPyReadCounter_files.baseName}' -c '$chromosome_files'
    """

}

    reads_ch      = channel.fromPath(params.reads, checkIfExists: true)
    barcodes_ch   = channel.fromPath(params.barcodes, checkIfExists: true)
    novoindex_ch  = channel.fromPath(params.novoindex, checkIfExists: true)
    gtf_ch        = channel.fromPath(params.gtf, checkIfExists: true)
    genome_ch     = channel.fromPath(params.genome, checkIfExists: true)

workflow {

    // run flexbar to trim adapter sequences
   
    flexbar_ch = FLEXBAR( reads_ch ) 

    //run pyBarcodeFilter with flexbar trimmed reads and barcode sequences

    demultiplex_ch = DEMULTIPLEX( flexbar_ch, barcodes_ch )

    //Collapse the demultiplexed files to remove duplicates using pyFastqDuplicateRemover, flatten the outputs of demultiplexing to stage the inputs as seven seperate files

    collapse_input_ch = demultiplex_ch.flatten()

    collapse_ch = COLLAPSE( collapse_input_ch )

    //align the processed reads to the pregenerated Novoindex using Novoalign

    novoalign_ch = NOVOALIGN( novoindex_ch, collapse_ch)

    //generation of hit tables with pyReadCounters.py from the aligned reads

    RUNPYREADCOUNTERS(gtf_ch, novoalign_ch)

    mapped_ch = RUNSECONDPYREADCOUNTERS(gtf_ch, novoalign_ch)

    //generation of chromosome files for usage by the subsequent steps

    chromosome_ch = CHROMOSOMELENGTH( genome_ch )

    //generation of a readable coverage file using pyGTF2sgr.py

    COVERAGE( chromosome_ch, mapped_ch )

    //generation of readable bedgraph file using pyGTF2bedgraph.py

    BEDGRAPH( chromosome_ch, mapped_ch )

}

workflow.onComplete = {
    //Print message on pipeline completion
    println "ProtoNF_novo.nf run Successful"
    println "Command run: $workflow.commandLine"
    println "Results published to pipeline directory under /results/"
}