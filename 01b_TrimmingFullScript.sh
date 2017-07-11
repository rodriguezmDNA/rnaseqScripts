#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Trim sequencing files in two steps.
## Three options are offered to test each step of the process.
## If there is adapter contamination (likely to be determined by the output of FastQC on raw files or using minion) \
## The option C can be used.
## If one is certain there is no adapter, option B works as well.
# Last modified
# 2016.06.28


# Set options:
# A. Remove only barcode (first 8bp) [fastx]
# B. Remove barcodes and trim based on quality, complexity (trims repeated bases and Ns) [fastx & reaper]
# C. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]
optionTrim=C



tabu='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
seqAdapt='ATCTCGTATGCCGTCTTCTGCTTG'

#tabu="ACACGTCTGAACTCCAGTCACACTCAGGTATCTCGTATGCCG,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,ATCTCGTATGCCGTCTTCTGCTTG,GAAGAAGAAGAAA,CGGAAGAGCACACGTCTGAACTCCAGTCAC,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,CACGTCTGAACTCCAGTCACTACCATTGATCTCGTATGCCGT,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#seqAdapt="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
## Possible contaminants

# ATCTCGTATGCCGTCTTCTGCTTG small RNA 3' Adapter
# GAAGAAGAAGAAA from minion. Run: minion search-adapter -i <file.fastq.gz>

# CGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter, extended


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
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).$optionTrim.log  
##
echo `touch $logPath` #Create file to fill with log data
echo 2>&1 | tee $logPath ## Empty the file contents each time the script is run 
echo "Started on `date`" 2>&1 | tee -a $logPath

######
echo "Trimming option $optionTrim" 2>&1 | tee $logPath

########################################################################################################################
########################################################################################################################

SeqDir="$DIR/RawData/SALK-012253-1-*.fastq.gz"
#SeqDir="$DIR/Data/demultiplexed/*.fq.gz"


## Remove first 8bp and sequences with low quality (< 35) with fastx_trimmer
## In this special case, the seq. facility removes adapter. Reaper will only trim reads based on 
#### quality and complexity.  


## Choose flags for fastx trimmer and reaper
if [ $optionTrim == "A" ]; then   ##Only trim
	echo 'A. Remove only barcode (first 8bp) [fastx]' 2>&1 | tee -a $logPath;
	ToDir="$DIR/01_trimmed/"  #ToDir="$DIR/01_trimmed/A_trimBC_noReap"
	fastxParams="-v -f 9 -Q33 -z -o $ToDir/$base.nobc.fq.clean.gz" #for fastx_trimmer, keep consistency with how reaper outputs clean files (.clean)

else 
	if [ $optionTrim == "B" ]; then ##Coupled with reaper
	echo 'B. Remove barcodes and trim based on quality, complexity (trims repeated bases and Ns) [fastx & reaper]' 2>&1 | tee -a $logPath;
	
	## Set params
	ToDir="$DIR/01_trimmed/" #ToDir="$DIR/01_trimmed/B_trimBC_noAdapt_Reap" 
	fastxParams="-v -f 9 -Q33" # If coupled with reaper
	reapParams='-geom no-bc -3pa "" -dust-suffix 6/ACTG -nnn-check 1/1 -qqq-check 35/10 -clean-length 30 -tri 20 -polya 6' #for reaper; removed --fastqx-out
else 
	if [ $optionTrim == "C" ]; then ##Coupled with reaper
	echo 'C. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]' 2>&1 | tee -a $logPath; 
	## Set params
	ToDir="$DIR/01_trimmed/" #ToDir="$DIR/01_trimmed/C_trimBC_trimAdapt_Reap" 
	fastxParams="-v -f 9 -Q33" # If coupled with reaper

	reapParams='-geom no-bc -dust-suffix 6/ACTG --noqc -nnn-check 1/1 -qqq-check 35/10 -clean-length 30 -tri 20 -polya 6' #for reaper; took away --fastqx-out
else 
	## Exit if no option selected
	echo 'Wrong option, exiting program' 2>&1 | tee -a $logPath;
	exit
fi; fi; fi; ## Close every if opened
########### Done setting up parameters


## Crete out dir
echo `mkdir -p $ToDir`


############
### Process sequencing files

for file in $SeqDir
do # do something on "$file"
 base=`basename $file`; # Remove dir paths
 base=${base%%.fastq.gz};  # Remove suffixes
 echo "Trimming $base" 2>&1 | tee -a $logPath # Print file being processed

 ##############################
 # First unzip to standar output
 ## Remove first 8bp and sequences with low quality (< 33) using fastx_trimmer
 ## Trim adapter, low complexity/quality sequences with reaper
 ## stdout and stderr to logfile while still display on screen
 

## Probably there's another fancier way to do this: 
### Run trimmer with/without reaper:
if [ $optionTrim == "A" ]; then   ##Only trim
	echo 'A. Remove only barcode (first 8bp) [fastx]' #2>&1 | tee -a $logPath;
	#ToDir="$DIR/02_trimmed.trimBC_noReap"
	#fastxParams="-v -f 9 -Q33 -z -o $ToDir/$base.nobc.fq.gz" #for fastx_trimmer
###
	echo `gunzip -dc $file |  fastx_trimmer $fastxParams` #2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test
else 
	if [ $optionTrim == "B" ]; then ##Coupled with reaper
	echo 'B. Remove barcodes and trim based on quality, complexity (trims repeated bases and Ns) [fastx & reaper]';
	## Set params
	#ToDir="$DIR/02_trimmed.trimBC_noAdapt_Reap" 
	#fastxParams="-v -f 9 -Q33" # If coupled with reaper
	#reapParams='-geom no-bc -3pa "" -dust-suffix 6/ACTG --fastqx-out -nnn-check 1/1 -qqq-check 35/10 -clean-length 15 -tri 20 -polya 6' #for reaper
###
	echo `gunzip -dc $file |  fastx_trimmer $fastxParams | reaper -basename $ToDir/$base $reapParams` #2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test

else 
	if [ $optionTrim == "C" ]; then ##Coupled with reaper
	echo 'C. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]' #2>&1 | tee -a $logPath; 
	echo 'Using' $tabu 'as tabu sequences'  2>&1 | tee -a $logPath;
	echo 'Using' $seqAdapt 'as adapter sequences' 2>&1 | tee -a $logPath;
	## Set params
	#ToDir="$DIR/02_trimmed.trimBC_trimAdapt_Reap" 
	#fastxParams="-v -f 9 -Q33" # If coupled with reaper
	#reapParams='-geom no-bc -3pa AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -dust-suffix 6/ACTG --fastqx-out -nnn-check 1/1 -qqq-check 35/10 -clean-length 15 -tri 20 -polya 6' #for reaper
####
	echo `gunzip -dc $file | fastx_trimmer $fastxParams | reaper -3pa ATCTCGTATGCCGTCTTCTGCTTG -tabu GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -basename $ToDir/$base $reapParams` #2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test
fi; fi; fi; ## Close every if opened
#### 


echo "----------------------------------------------------------" #2>&1 | tee -a $logPath
done >&2 | tee -a  $logPath ## Only errors!



## Old arguments, if having adapter contamination.
## Previous to 16.04.18 version. Used possible adapter contaminants for BRAD-Seq
#### Trim adapter, low complexity/quality sequences with reaper
# `gunzip -dc $file | ../Software/fastx/fastx_trimmer -f 9 -v -Q 33 \
#                   | ../Software/Kraken/src/reaper -basename $ToDir/$base \
#	-geom no-bc -3pa GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC,GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG \
#        -tabu GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#        -mr-tabu 14/2/1 -3p-global 6/1/0/0 -3p-prefix 11/2/1 -3p-head-to-tail 1 \
#        -dust-suffix 7/ACTG --fastqx-out -nnn-check 1/1 -qqq-check 35/10 -clean-length 15 -polya 10`

########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $fastxParams" 2>&1 | tee -a $logPath
echo -e "\nParameters used: $reapParams" 2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath

