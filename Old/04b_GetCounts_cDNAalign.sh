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


## Read
inputFolder="bwt1_cDNA" #Or: bwt1_genome,bwt2_genome,bwt2_cDNA


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
echo -e "Getting gene counts from a cDNA aligned files \n" 2>&1 | tee -a $logPath

#######################################################################################################################

#inputFolder="$DIR/03b_deduplicated/"
#
outDir=$DIR/04_Counts/$inputFolder
countsDir=$outDir

echo `mkdir -p $countsDir` 2>&1 | tee -a $logPath;

## Run htseq 
for each in `ls $DIR/03_alignment/$inputFolder`; do
	Name=$(basename $each .bam)
	bamFile=$DIR/03_alignment/$inputFolder/$each 
	echo "Processing $each";
	#echo $bamFile
	# The next step will:
	
	# Convert BAM to SAM | limit to the 3rd column. Sort and count the number of occurences of a gene. Swap the columns. Save to file.
	
	echo `samtools view -q 30 $bamFile | cut -f 3 | sort | uniq -c | awk '{ t=$1 ; $1=$2; $2=t; print }'  > $countsDir/$Name.counts.txt`;
done 2>&1 | tee -a $logPath


## Count stats:
#logCount=$DIR/logs/Counts #Create log folder if it doesn't exist
#echo `mkdir -p $logCount`;

#for each in `ls $countsDir`; do
#        countName=$(basename $each .txt)
#        countFile=$countsDir/$each
#	echo `awk '/^__/ { print $0 }' $countFile  > $logCount/$countName.stat.log` #Get all lines starting with double underscore into a new file
#    echo `awk '!/^__/ { print $0 }' $countFile > $outDir/$countName.txt` #Counts file without stats
#done 2>&1 | tee -a $logPath


   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done
