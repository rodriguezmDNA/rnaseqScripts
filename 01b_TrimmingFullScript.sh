#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Trim sequencing files in two steps.
# If there is adapter contamination (likely to be determined by the output of FastQC on raw files or using minion) 


## If one is certain there is no adapter, leave the tabu and seqAdapt blank
# ie: tabu='';seqAdapt=''

## If one is certain the barcodes have been removed (first 6-8 bp) set the -f parameter for fastx to 1.

#####
# Last modified
# 2017.07.jrm
####################


#####
# Set options:
# A. Remove only barcode (first 8bp) [fastx]
# B. Remove barcodes. Remove adapters. Trim based on quality, complexity (trims repeated bases and Ns) [fastx & reaper]

optionTrim=B

# Barcode trim options.
ntStart=9 #Non-zero positive integers
## eg: A value of 9 removes the first 8 nt. 
## In some cases the sequencing facility removes barcodes, if so, set ntStart=1


# Reaper options
tabu='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' #Reaper recommends to use the 5' adapter sequence as tabu. If more than one contaminant is suspected, they can be inserted as a comma separated list: tabu='AAA,GGG'
seqAdapt='ATCTCGTATGCCGTCTTCTGCTTG' # 3' adapter. See below.
####################


#tabu="ACACGTCTGAACTCCAGTCACACTCAGGTATCTCGTATGCCG,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,ATCTCGTATGCCGTCTTCTGCTTG,GAAGAAGAAGAAA,CGGAAGAGCACACGTCTGAACTCCAGTCAC,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,CACGTCTGAACTCCAGTCACTACCATTGATCTCGTATGCCGT,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#seqAdapt="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
## Possible contaminants

# ATCTCGTATGCCGTCTTCTGCTTG small RNA 3' Adapter
# GAAGAAGAAGAAA from minion. Run: minion search-adapter -i <file.fastq.gz>

# CGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter, extended

#####
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

SeqDir="$DIR/RawData/*.fastq.gz"

## Remove first 8bp and sequences with low quality (< 33) with fastx_trimmer



## Choose flags for fastx trimmer and reaper
if [ $optionTrim == "A" ]; then   # Only trim
	echo 'A. Remove only barcode (first 8bp) [fastx]' 2>&1 | tee -a $logPath; #Describe

	# Set parameters
	ToDir="$DIR/01_trimmed/"  #ToDir="$DIR/01_trimmed/A_trimBC_noReap"
	fastxParams='-v -Q33 -z' #for fastx_trimmer, keep consistency with how reaper outputs clean files (.clean). This will help on subsequent steps.
else 
	if [ $optionTrim == "B" ]; then # Coupled with reaper
	echo 'C. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]' 2>&1 | tee -a $logPath; 
	## Set params
	ToDir="$DIR/01_trimmed/" #ToDir="$DIR/01_trimmed/C_trimBC_trimAdapt_Reap" 
	fastxParams='-v -Q33' # If coupled with reaper
	reapParams='-geom no-bc -dust-suffix-late 10/ACTG -dust-suffix 10/ACTG --noqc -nnn-check 1/1 -qqq-check 33/10 -clean-length 30 -tri 20 -polya 5 --bcq-late' #for reaper; took away --fastqx-out
else 
	## Exit if no option selected
	echo 'Wrong option, exiting program' 2>&1 | tee -a $logPath;
	exit
fi; fi; ## Close every if opened
########### Done setting up parameters


## Crete out dir
echo `mkdir -p $ToDir`


############
### Process sequencing files

# For each file in the selected folder
## Fix the names
## Depending on the option selected:
## A) Just trim the first n letters.
## B) Trim and remove adapters. Clean sequences based on quality/complexity. 
## Write the output to a file in a dedicated folder. 
## Errors and screen outputs go to both screen and file.

for file in $SeqDir
do # do something on "$file"
 base=`basename $file`; # Remove dir paths
 base=${base%%.fastq.gz};  # Remove suffixes
 echo "Trimming $base" 2>&1 | tee -a $logPath # Print file being processed

 ##############################
 # First unzip to standar output
 ## Remove first n-bp and sequences with low quality (< 33) using fastx_trimmer
 ## Trim adapter, low complexity/quality sequences with reaper
 ## stdout and stderr to logfile while still display on screen
 
## Probably there's another fancier way to do this: 
### Run trimmer with/without reaper:
if [ $optionTrim == "A" ]; then   ##Only trim
	echo 'A. Remove only barcode (first 8bp) [fastx]' #2>&1 | tee -a $logPath;
	
	## Unzip to stdout | Trim the first n-bases and save output  
	gunzip -dc $file | fastx_trimmer -f $ntStart -Q33 $fastxParams -o $ToDir/$base.nobc.fq.clean.gz 2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test
else 
	if [ $optionTrim == "B" ]; then ##Coupled with reaper
	echo 'B. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]' #2>&1 | tee -a $logPath; 
	echo 'Using' $tabu 'as tabu sequences'  2>&1 | tee -a $logPath;
	echo 'Using' $seqAdapt 'as adapter sequences' 2>&1 | tee -a $logPath;
		
	####
	# Unzip the file, cut the first n bp, trim adapters and low quality segments.
	gunzip -dc $file |  fastx_trimmer -f $ntStart $fastxParams | reaper -3pa $seqAdapt -tabu $tabu -basename $ToDir/$base $reapParams 2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test
fi; fi; ## Close every if opened
#### 


echo -e "\n ----------------------------------------------------------" 2>&1 | tee -a $logPath
done >&2 | tee -a  $logPath ## Only errors!


########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $fastxParams -f $ntStart" 2>&1 | tee -a $logPath
echo -e "\nParameters used: $reapParams -3pa $seqAdapt -tabu $tabu" 2>&1 | tee -a $logPath
echo -e "\nExecution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\nDone `date`"  2>&1 | tee -a $logPath

