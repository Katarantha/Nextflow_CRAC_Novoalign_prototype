# CRAC-pipeline-Novoalign-Prototype

Author: Connor Hudson
Date Finalised: 8th March 2023

Introduction
This script utilised the pyCRAC package developed by Sander Granneman and Novoalign software to take CRAC RNA protein interaction data and process, align and subsequently process into bedgraph and sgr files.

Dependancies
This pipeline was developed on dependancy verions:
Nexflow: version 21.10.6

Flexbar: version 3.5.0

PyCRAC software tools:
pyBarcodeFilter.py             : version 3.1
pyFastqDuplicateRemover.py     : version 0.0.5
pyReadCounters.py              : version 0.5.5
pyCalculateChromosomeLengths.py: version 0.0.5
pyGTF2sgr.py                   : version 0.1.0
pyGTF2bedGraph.py              : version 0.1.0


Novoalign: version 2.07.00
note this software is proprietary and does require a payed license

Usage:
This pipeline requires 5 user inputs:
Sample Data        : specified by --reads on command line
Barcodes           : specified by --barcodes 
Novoalign Index    : specified by --novoindex
GTF Annotation file: specified by --gtf
Genome Fasta file  : specified by --genome

Outputs
This pipelien will out 5 types of relevant outputs to the directory the pipeline is run in under a new directory named results/

results/Alignments: outputs from the Novoalign process

results/hitTable: outputs from the first mapping of genomic features process

results/hitTable2: outputs from the second mapping of genomic features process

results/bedgraph: outputs from generation of bedgraph files from mapped results

results/sgrCoverage: results from generation of sgr Coverage files from the mapped results



