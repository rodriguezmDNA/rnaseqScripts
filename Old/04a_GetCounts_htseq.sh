#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Get counts from ht-seq
# Last modified
# jrm_07_2017

## Only for genome-aligned files.
## Run the associated R file after this to generate a matrix of counts. This can be then used for DGE with your favourite program (DESeq, edgeR, limma, poissonSeq).


# Before starting, check that the annotation files and the fasta files with the genome have matching names.
## In the command line manually change the names for plastid genomes in the Arabidopsis thaliana annotation (TAIR v10):
# in vim, type ":" to enter command mode and pase:
# %s/ChrM/mitochondria/g
# %s/ChrC/chloroplast/g
####################

# Get counts from BAM files (output from tophat) using HTSeq-count
## Record time
start_time=`date +%s`

## Get absolute path
#Get the full path of the current script: Go back one folder (from scripts to main)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

## Keep track of errors and outputs in a log.
logDir=$DIR/logs #Create log folder if it doesn't exist
if [ ! -d $logDir ]; then echo `mkdir -p $logDir`; fi
## 
######
##Uses the same name of the script but adds the log extension
## Two ways: split string by dots, keep the first part
#logPath=$DIR/logs/$(basename $BASH_SOURCE | cut -d . -f1).log 
# Removes specific extension:
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).log

##
echo `touch $logPath` #Create file to fill with log data
echo "Started on `date`" 2>&1 | tee $logPath
######

######
echo -e "Getting gene counts with HTSeq-count  \n" 2>&1 | tee -a $logPath

#######################################################################################################################

#Index suffix
annotName='TAIR10_withTransposons' #Same as BuildIndex script
annotExt='.gff' #Can be gff or gtf
#First, recreate fasta files from indices
annotationPath=$DIR"/meta/$annotName$annotExt"

## Read
inputFolder="bwt1_genome"
#inputFolder="$DIR/03b_deduplicated/"
#
outDir=$DIR/04_Counts/$inputFolder
countsDir=$outDir/htseqraw
echo `mkdir -p $countsDir` 2>&1 | tee -a $logPath;

## Set parameters for tophat
# cut -f3 TAIR10_withTransposons.gff | sort | uniq # To see the optios for  --type:
biotype="gene" #mRNA, CDS, 
attribute="ID" #ID, Parent
mode="intersection-strict" #"union"
htseqParams="-f sam --stranded yes --mode $mode --idattr $attribute --type $biotype" #--format=sam

## Run htseq 
for each in `ls $DIR/03_alignment/$inputFolder`; do
	Name=$(basename $each .bam)
	bamFile=$DIR/03_alignment/$inputFolder/$each 
	echo "Processing $each";
	# The next step will:
	# Convert BAM to SAM | sort the aligned file in SAM format | use htseq-count to get the number of reads per biotype (gene,exon,etc).
	###### Test the first 15 alignemnts
	#echo `samtools view -h $bamFile | samtools sort -O sam | head -n 15| htseq-count $htseqParams - $annotationPath > $countsDir/$Name.counts.txt`;
	# Run the whole thing:
	echo `samtools view -h $bamFile | samtools sort -O sam |  htseq-count $htseqParams - $annotationPath > $countsDir/$Name.counts.txt`;
done 2>&1 | tee -a $logPath


## Count stats:
logCount=$DIR/logs/Counts #Create log folder if it doesn't exist
echo `mkdir -p $logCount`;

for each in `ls $countsDir`; do
        countName=$(basename $each .txt)
        countFile=$countsDir/$each
	echo `awk '/^__/ { print $0 }' $countFile  > $logCount/$countName.stat.log` #Get all lines starting with double underscore into a new file
    echo `awk '!/^__/ { print $0 }' $countFile > $outDir/$countName.txt` #Counts file without stats
done 2>&1 | tee -a $logPath


   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: "  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done
