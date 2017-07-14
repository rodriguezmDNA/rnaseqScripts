#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Create index file for bowtie2 mapping
# Last modified
# 2016.06.18s


# Why another script?
##
## Bowtie 2's command-line arguments and genome index format are both different 
#  from Bowtie 1's.
## When to use bowtie1 or bowtie2?
# If reads are >50bp, use bowtie2
# If reads are shorter than 50bp, use bowtie1 (aka, bowtie)

# Don't use masked genomes

# Set options:
# bwt1 -> To use bowtie1 
# bwt2 -> To use bowtie1
optionBowtie='bwt1'


## Record time
start_time=`date +%s`

## Get absolute path
#Get the full path of the current script: Go back one folder (from scripts to main)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

## Keep track of errors and outputs in a log.
logDir=$DIR/logs #Create log folder if it doesn't exist
if [ ! -d $logDir ]; then echo `mkdir -p $logDir`; fi

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
echo -e "Index Creation with option $optionBowtie \n" 2>&1 | tee $logPath


#### Functions
checkDir () {
	if [ ! -d $toDir ]; then     
		echo "Creating directory $toDir" 2>&1 | tee -a $logPath;
		mkdir -p $toDir;
	else  
		echo "Directory $toDir found" 2>&1 | tee -a $logPath;
	fi;}

           
checkIndex  () {
	if [ ! -f $toDir/$pattern ]; then
		echo -e "Building index\n" 2>&1 | tee -a $logPath;
		#echo $buildIndex
		## Build index:
		$bwt $buildIndex 2>&1 | tee -a $logPath;
	else
		echo -e "An index file was found, not building new index \n" 2>&1 | tee -a $logPath;
	fi;}

########################################################################################################################
########################################################################################################################

#In params:
genomeFasta=$DIR'/meta/genome_ath_TAIRv10.fa'
#genomeFasta='300.fa' ## Only for testing
cdnaFasta=$DIR'/meta/cdna_ath_TAIRv10.fa'
suffixOut='athIndex'

pattern="ath"
# Out params:
genomeToPath=$DIR"/"$pattern"Genome/"$optionBowtie"_genome"
cdnaToPath=$DIR"/"$pattern"Genome/"$optionBowtie"_cdna"



## Choose flags for Bowtie or Bowtie 2
if [ $optionBowtie == "bwt1" ]; then   
	echo "Bowtie1"
	# Set Global Parameters
	pattern=*.ebwt
	bwt=bowtie-build
	
	## Genome
	echo -e "\n Genome Index"
	toDir=$genomeToPath
	checkDir $toDir
	# Check if index exists
	buildIndex="$genomeFasta $toDir/$suffixOut"
	checkIndex $toDir $pattern $buildIndex $bwt

	# cDNA
	echo -e "\n Transcriptome Index"
	toDir=$cdnaToPath
	checkDir $toDir
	# Check if index exists
	buildIndex="$cdnaFasta $toDir/$suffixOut"
	checkIndex $toDir $pattern $buildIndex $bwt

else 
## Bowtie 2
if [ $optionBowtie == "bwt2" ]; then
	echo "Bowtie2"
	# Set Global Parameters
	pattern=*.bt2
	bwt=bowtie2-build

	## Genome
	echo -e "\n Genome Index"
	toDir=$genomeToPath
	checkDir $toDir
	# Check if index exists
	buildIndex="$genomeFasta $toDir/$suffixOut"
	checkIndex $toDir $pattern $buildIndex $bwt


	# cDNA
	echo -e "\n Transcriptome Index"
	toDir=$cdnaToPath
	checkDir $toDir
	# Check if index exists
	buildIndex="$cdnaFasta $toDir/$suffixOut"
	checkIndex $toDir $pattern $buildIndex $bwt
	



	
fi;fi;

# Set parameters
#

## Check if directory and/or index exist

# Directory


## Testing for index


   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $optionBowtie $genomeFasta $cdnaFasta"  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done
