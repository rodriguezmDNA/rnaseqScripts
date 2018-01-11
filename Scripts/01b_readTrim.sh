#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Barcode trimming and cleaning.
# Last modified
# 2017.10.jrm

#### Some predefined options. Modify at your own risk
display_help() {
	echo "Wrapper to trim and clean RNASeq data." >&2
	#echo "If the script detects an existing index it won't create a new one." >&2
	#echo "." >&2
	echo ""
	echo "Usage: $0 [options...] " >&2
    echo ""
    echo
    echo "   -h            this useful help"
    echo "   -m            mode: A) Barcode removal (First n bp)"
    echo "                       B) Barcode & trim (Uses reaper)"
    echo "   -n           length of barcode (removes first n bases)."
    echo "                ie -n 9 will remove the first 8bp."
    echo "                Default value: 1; Use 8 for processing BRAD-Seq data"
    echo "   -p           Path to a file with parameters for reaper - See reaper documentation "
    echo "   -t           tabu sequence - See reaper documentation "
    echo "   -a           adapter sequence - See reaper documentation"
    echo "                -t '' or -a '' sets sequence as empty otherwise uses default"
    echo "                -t 'AAA,GGG' for multiple sequences"
    echo "Assumptions:"
    echo "* A directory with Raw reads (with extension \".fq or .fastq\" on a folder named 'RawData/')"
    echo "* Uses fastx_trimmer and/or reaper. They must be accessible through the PATH')"
    echo ""
	echo "---	Default values:"
	echo "fastx_trimmer:"
	echo "-n: 1"
	echo "reaper:"
    echo "parameters: -geom no-bc -dust-suffix-late 10/ACTG -dust-suffix 10/ACTG --noqc -nnn-check 1/1 -qqq-check 33/10 -clean-length 30 -tri 20 -polya 5 --bcq-late"
    echo "tabu [-t]:  GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 
    echo "adapter [-a]: ATCTCGTATGCCGTCTTCTGCTTG"
    echo ""
    echo -e "\n don't panic! "
    echo " github: rodriguezmDNA "
    echo " Last modified: 2017.10 - jrm - @ UC Davis "
    
    
    exit 1
}

if [ $# -eq 0 ]
  then
    display_help
fi

#####
# Set options:
# A. Remove only barcode (first 8bp) [fastx]
# B. Remove barcodes. Remove adapters. Trim based on quality, complexity (trims repeated bases and Ns) [fastx & reaper]


###### fastx_trim options
####################
ntStart=1 #Non-zero positive integers

####################

###### Reaper options
####################
tabu='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' #Reaper recommends to use the 5' adapter sequence as tabu. 
#If more than one contaminant is suspected, they can be inserted as a comma separated list: tabu='AAA,GGG'
seqAdapt='ATCTCGTATGCCGTCTTCTGCTTG' # 3' adapter. See below.


#tabu="ACACGTCTGAACTCCAGTCACACTCAGGTATCTCGTATGCCG,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,ATCTCGTATGCCGTCTTCTGCTTG,GAAGAAGAAGAAA,CGGAAGAGCACACGTCTGAACTCCAGTCAC,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,CACGTCTGAACTCCAGTCACTACCATTGATCTCGTATGCCGT,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#seqAdapt="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
## Possible contaminants

# ATCTCGTATGCCGTCTTCTGCTTG small RNA 3' Adapter
# GAAGAAGAAGAAA from minion. Run: minion search-adapter -i <file.fastq.gz>

# CGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC #TruSeq adapter, extended

while getopts ':hm:n:t:a:p:' option; do
  case "$option" in
    h) display_help
       exit
       ;;
    m) optionTrim=$OPTARG
		if [ $optionTrim == 'A' ] || [ $optionTrim == 'B' ]
		    then
		    echo ""
		else
		    echo "Mode not valid, use either -m A or -m B" >&2
		    display_help
		fi
		;;
    t) tabu=$OPTARG
       ;;
    a) seqAdapt=$OPTARG
       ;;
    n) ntStart=$OPTARG
       ;;
    p) reapParams=$OPTARG
	   ;;
    :) printf "missing value for -%s\n" "$OPTARG" >&2
       display_help
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       display_help
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))


######### Validate parameters
############################################

if (( $ntStart < 1 ))
	then
	echo "-n $ntStart unvalid. Please use a positive, non-zero value. Default: 1"
	display_help
fi


### If empty, set reaper to default parameters
if [[  -z  $reapParams ]]; then
	echo "Reaper parameters are empty, using default parameters" >&2

	reapParams='-geom no-bc -dust-suffix-late 10/ACTG -dust-suffix 10/ACTG --noqc -nnn-check 1/1 -qqq-check 33/10 -clean-length 30 -tri 20 -polya 5 --bcq-late' #for reaper; took away --fastqx-out    
 
else
	### Read parameters
	echo "Reading $reapParams"
	reapParams=`cat $reapParams`
fi

#echo $reapParams
#echo "`reaper $paramsFile`"


# ### If empty, set to default parameters
# if [[  -z  $paramsFile ]] ; then
# 	echo "Bowtie parameters are empty, using default parameters" >&2
# 	case $optionBowtie in
    
#     bwt1)
# 		# Bowtie 1 parameters I use
# 		bwtParams="-a --best --strata -n 1 -m 1 -p 4 --sam --tryhard"
#         ;;
#     bwt2)
# 		bwtParams="--time -p3 -k 10 --sensitive-local --no-unal"
# 		;;
#     # *)
#     #     echo "Option not valid" >&2
#     #     display_help
#     #     ;;
# 	esac
# else
# 	### Read parameters
# 	bwtParams=`cat $paramsFile`
# fi


############################################
######### Start of the Analysis
############################################


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
#logBow=$DIR/logs/bowStats"_"$optionAlign$optionBowtie #Create log folder if it doesn't exist
#if [ ! -d $logBow ]; then echo `mkdir -p $logBow`; fi
######

## Use the same script's name but add the log extension
## Two ways: split string by dots, keep the first part
# logPath=$DIR/logs/$(basename $BASH_SOURCE | cut -d . -f1).log # What if the script name has more dots?
## Removes specific extension:
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh)_option$optionTrim.log

## Create a log file and start writing.
echo `touch $logPath` #Create file to fill with log data
echo "Started on `date`" 2>&1 | tee $logPath


########################################################################################################################
########################################################################################################################

echo "Trimming option $optionTrim" 2>&1 | tee $logPath

SeqDir="$DIR/RawData/*.f*.gz"


## Choose flags for fastx trimmer and reaper
if [ $optionTrim == "A" ]; then   # Only trim
	echo 'A. Remove only barcode (first 8bp) [fastx]' 2>&1 | tee -a $logPath; #Describe

	# Set parameters
	ToDir="$DIR/01_trimmed/"  #ToDir="$DIR/01_trimmed/A_trimBC_noReap"
	fastxParams='-v -Q33 -z' #for fastx_trimmer, keep consistency with how reaper outputs clean files (.clean). This will help on subsequent steps.
else 
	if [ $optionTrim == "B" ]; then # Coupled with reaper
	echo 'B. Remove barcodes and adapters. Trim reads based on quality and complexity [fastx & reaper]' 2>&1 | tee -a $logPath; 
	## Set params
	ToDir="$DIR/01_trimmed/" #ToDir="$DIR/01_trimmed/C_trimBC_trimAdapt_Reap" 
	fastxParams='-v -Q33' # If coupled with reaper
	
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
	gunzip -dc $file |  fastx_trimmer -f $ntStart -Q33 $fastxParams | reaper -3pa $seqAdapt -tabu $tabu -basename $ToDir/$base $reapParams 2>&1 | tee -a $logPath # Use head -n 4000 in the middle for a quick test
fi; fi; ## Close every if opened
#### 


echo -e "\n ----------------------------------------------------------" 2>&1 | tee -a $logPath
done >&2 | tee -a  $logPath ## Only errors!



   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`
# reaper -3pa $seqAdapt -tabu $tabu $reapParams
echo -e "\n\t Parameters used: \n
mode: $optionTrim (if mode A, ignore reaper arguments)
fastx_trimmer -f $ntStart -Q33 $fastxParams
reaper: -3pa $seqAdapt -tabu $tabu $reapParams"  2>&1 | tee -a $logPath
echo -e "\n\texecution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n\tDone `date`"  2>&1 | tee -a $logPath
##Done
