#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Create index file for bowtie2 mapping
# Last modified
# 2017.07.jrm

# !!!!!!!!!!!!!!!!!!!!!!!!!!
# Don't use masked genomes


## Bowtie 2's command-line arguments and genome index format are both different from Bowtie 1's.
## When to use bowtie1 or bowtie2?
# If reads are >50bp, use bowtie2
# If reads are shorter than 50bp, use bowtie1 (aka, bowtie)

############
# Set options:
# {bwt1/bwt2} -> Align with bowtie1 OR bowtie2
# {cdna/genome} -> Align to cDNA or genome

optionBowtie='bwt1' #bwt1 or bwt2
optionAlign='genome' #genome or cDNA
################################################################################################


############
## Record time
start_time=`date +%s`

## Get absolute path
#Get the full path of the current script: Go back one folder (from scripts to main)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

## Keep track of errors and outputs in a log.
logDir=$DIR/logs #Create log folder if it doesn't exist
if [ ! -d $logDir ]; then echo `mkdir -p $logDir`; fi
## 
logBow=$DIR/logs/bowStats"_"$optionAlign$optionBowtie #Create log folder if it doesn't exist
if [ ! -d $logBow ]; then echo `mkdir -p $logBow`; fi
######
##Uses the same name of the script but adds the log extension
## Two ways: split string by dots, keep the first part
#logPath=$DIR/logs/$(basename $BASH_SOURCE | cut -d . -f1).log 
# Removes specific extension:
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).$optionAlign$optionBowtie.log

##
echo `touch $logPath` #Create file to fill with log data
echo "Started on `date`" 2>&1 | tee $logPath
################################################################################################

######
echo -e "Align to $optionAlign with $optionBowtie \n" 2>&1 | tee -a $logPath


#### Functions
bow1 () {
	echo -e 'Align with bowtie 1 \n' 2>&1 | tee -a $logPath
	
	for fqFile in `ls $inputFolder/*$inSuffix`; do
		name=`echo $fqFile | perl -pe "s/.+\/|\.lane.clean.gz$//g"`; ##Remove suffix from file name
		echo $name;
		
		#For tests use | head -n 400 | in between gzip and bowtie		

		# Run in two steps to be able to display and save the output. 
		# unzip the file, use bowtie1 with the parameters chosen. Save as SAM. 
		gunzip -dc $fqFile | bowtie $bwtParams $idxFile - -S $outPath/$name.sam 2>&1 | tee $logBow/$name.log;
		# read the SAM file convert to BAM.
		samtools $toBAM $outPath/$name.sam -o $outPath/$name.bam 2>&1 | tee -a $logPath;
		# Delete the SAM file
		rm $outPath/$name.sam;
	done;

	}

           
bow2  () {
	echo -e 'Align with bowtie 2 \n' 2>&1 | tee -a $logPath
	
	for fqFile in `ls $inputFolder/*$inSuffix`; do
		name=`echo $fqFile | perl -pe "s/.+\/|\.$inSuffix$//g"`; ##Remove suffix from file name
		echo $name;
		echo $bwtParams;
		echo $idxFile;
		echo $toBAM;

		#For tests use | head -n 400 | in between gzip and bowtie		

		gunzip -dc $fqFile | bowtie2 $bwtParams -x $idxFile - -S $outPath/$name.sam 2>&1 | tee $logBow/$name.log;
		samtools $toBAM $outPath/$name.sam -o $outPath/$name.bam - 2>&1 | tee -a $logPath; #For tests use | head -n 400 | in between gzip and bowtie		
		# Delete the SAM file
		rm $outPath/$name.sam;
	done; #>&2 | tee $logBow/$name.log;

	}

########################################################################################################################
########################################################################################################################

#Index suffix
suffixOut='athIndex' #Same as BuildIndex script


## Read
inSuffix="*clean*"
inputFolder="$DIR/01_trimmed"


## Use Out params from previous script for index file
outPath=$DIR/"03_alignment"/$optionBowtie"_"$optionAlign
echo `mkdir -p $outPath`


## cDNA or transcriptome
if [ $optionAlign == "genome" ]; then
	echo "Genome"  2>&1 | tee -a $logPath
	idxFile=$DIR/"athGenome"/$optionBowtie"_genome"/$suffixOut
	strandAlign=""

else
	echo "cDNA" 2>&1 | tee -a $logPath
	idxFile=$DIR/"athGenome"/$optionBowtie"_cdna"/$suffixOut
	strandAlign="--norc"
fi;


## Align
if [ $optionBowtie == "bwt1" ]; then   
	#bwt1Params="--time -p3 -k 10 --sensitive-local --no-unal"
	bwtParams="-a --best --strata -n 1 -m 1 -p 4 --sam --tryhard"
	# Samtools parameters
	toBAM="view -q 30 -b -S" 
	#Report reads with a mapq of at least 30
	bow1 $bwtParams $strandAlign $toBAM
else 

## Bowtie 2
if [ $optionBowtie == "bwt2" ]; then ##Coupled with reaper
	bwtParams="--time -p3 -k 10 --sensitive-local --no-unal"
	toBAM="view -h -q30 -b -S"
	bow2 $bwtParams $toBAM

fi;fi;


   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $bwtParams $strandAlign"  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done
