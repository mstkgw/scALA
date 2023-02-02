# scALA

This scALA pipeline is available for Linux.

###### Requirment ######
python >= 3.6
biopython
java >= 1.8
php = 7.x
NanoFilt >= 2.6.0
canu >= 1.9
flye >= 12.7
minimap2 >= 2.17 
samtools >= 1.9
blastn >= 2.9.0
jvarkit
MultiCSAR

######   Usage   ######
0. Write path of the required tools at the bigging of main.sh
1. Run main.sh as follows

    bash scALA.sh -i LongRead.fastq -o OutDir -g 3.1m -d 50

    options:
     -i, required. Input LongRead Fastq File (Not Compressed)
     -o, required. Output Directory
     -d, required. Threshold of Depth in Debiasing Step
     -g, required. Estimated Genome Size of Target Bacteria
     -r, optional. If you have a known reference genome.
     -s, optional. If you want to use short-read contig in the scaffolding step.

2. Get "*_scaffold.fasta" in output directory, and polish the sequence in the way you prefer.
