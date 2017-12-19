#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# Bowtie mapping 
# Last modified
# 2017.10.jrm

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

display_help() {
	echo "" >&2
	echo "Usage: $0 [option...] " >&2
    echo "Wrapper to run bowtie{1/2} on trimmed sequenced reads.  " >&2
    echo
    echo "   -h            this useful help"
    echo "   -v            bowtie version [bwt1/bwt2]"
    echo "   -r            reference [genome/cDNA]"
    echo "Optional: "
    echo "   -a            path to a text file with bowtie arguments"
    echo "   -s            suffix of the reference (default: \"RefIdx\")"
    echo
    echo "example: Scripts/makeIdx.sh -v bwt1 -r cDNA -s athIndex -a meta/myarguments.txt"
    echo
    echo "Assumptions:"
    echo "* A directory with trimmed reads (with extension \"clean\" on a folder named '01_trimmed/')"
    echo "* The index of the reference is in a 'RefIdx/' directory with format: "
	echo "  'bwt{v}_reference'. ie, that the directory 'RefIdx/bwt1_cdna/' exists and contains index files"
	echo "  built with the appropiate bowtie version"
    
    echo -e "\n don't panic! "
    echo " github: rodriguezmDNA "
    echo " Last modified: 2017.10 - jrm - @ UC Davis "
    
    
    exit 1
}

if [ $# -eq 0 ]
  then
    display_help
fi

#### Some predefined options. Modify at your own risk
## Read data. This assumes you have the 01_trimmed folder after trimming with reaper.

####

suffixOut='RefIdx' #Same as BuildIndex script
while getopts ':hv:r:a:s:' option; do
  case "$option" in
    h) display_help
       exit
       ;;
    v) optionBowtie=$OPTARG
		if [ $optionBowtie == 'bwt1' ] || [ $optionBowtie == 'bwt2' ]
		    then
		    echo ""
		else
		    echo "Bowtie version not valid" >&2
		    display_help
		fi
		;;
    r) optionAlign=$OPTARG
       ;;
    a) paramsFile=$OPTARG
       ;;
    s) suffixOut=$OPTARG
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

## Reference Option
if [ $optionAlign == 'cDNA' ] || [ $optionAlign == 'genome' ]
    then
    echo ""
else
    display_help
fi

### If empty, set to default parameters
if [[  -z  $paramsFile ]] ; then
	echo "Bowtie parameters are empty, using default parameters" >&2
	case $optionBowtie in
    
    bwt1)
		# Bowtie 1 parameters I use
		bwtParams="-a --best --strata -n 1 -m 1 -p 4 --sam --tryhard"
        ;;
    bwt2)
		bwtParams="--time -p3 -k 10 --sensitive-local --no-unal"
		;;
    # *)
    #     echo "Option not valid" >&2
    #     display_help
    #     ;;
	esac
else
	### Read parameters
	bwtParams=`cat $paramsFile`
fi


############################################
######### Start of the Analysis
############################################

echo "bwParams" $bwtParams

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
logBow=$DIR/logs/bowStats"_"$optionAlign$optionBowtie_$suffixOut #Create log folder if it doesn't exist
if [ ! -d $logBow ]; then echo `mkdir -p $logBow`; fi
######

## Use the same script's name but add the log extension
## Two ways: split string by dots, keep the first part
# logPath=$DIR/logs/$(basename $BASH_SOURCE | cut -d . -f1).log # What if the script name has more dots?
## Removes specific extension:
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).$suffixOut$optionAlign$optionBowtie.log

## Create a log file and start writing.
echo `touch $logPath` #Create file to fill with log data
echo "Started on `date`" 2>&1 | tee $logPath
############

inSuffix="*clean*" #Extension for clean reads (default from reaper)
inputFolder="$DIR/01_trimmed"


######### Functions
############################################

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
		name=`echo $fqFile | perl -pe "s/.+\/|\.lane.clean.gz//g"`; ##Remove suffix from file name
		echo $name;
		echo $bwtParams;
		echo $idxFile;
		echo $toBAM;

		#For tests use | head -n 400 | in between gzip and bowtie		

		gunzip -dc $fqFile | bowtie2 $bwtParams -x $idxFile - -S $outPath/$name.sam 2>&1 | tee $logBow/$name.log;
		#samtools $toBAM $outPath/$name.sam -o $outPath/$name.bam - 2>&1 | tee -a $logPath; #For tests use | head -n 400 | in between gzip and bowtie		
		## Get uniquely mapped reads:
		samtools view -h -q 30 $outPath/$name.sam | grep -v "XS:i:" | samtools sort - -o $outPath/$name.bam 2>&1 | tee -a $logPath; #For tests use | head -n 400 | in between gzip and bowtie		
		# Delete the SAM file
		rm $outPath/$name.sam*;
	done; #>&2 | tee $logBow/$name.log;

	}

######### Start of the Script
############################################


## Use Out params from previous script for index file
outPath=$DIR/"03_alignment"/$optionBowtie"_"$optionAlign_"$suffixOut"
echo `mkdir -p $outPath`


## cDNA or transcriptome
if [ $optionAlign == "genome" ]; then
	echo "Genome"  2>&1 | tee -a $logPath
	idxFile=$DIR/"RefIdx"/$optionBowtie"_genome"/$suffixOut
	strandAlign=""  
else
	echo "cDNA" 2>&1 | tee -a $logPath
	idxFile=$DIR/"RefIdx"/$optionBowtie"_cdna"/$suffixOut
	strandAlign="--norc"
fi;


## Align
if [ $optionBowtie == "bwt1" ]; then   
	
	# Some other options.
	#bwt1Params="--time -p3 -k 10 --sensitive-local --no-unal"
	#bwtParams="-a --best --strata -n 1 -m 1 -p 4 --sam --tryhard"
	
	## Samtools parameters
	#Report reads with a mapq of at least 30:
	toBAM="view -q 30 -b -S" 

	################
	# Run function
	echo "running bow1"
	bow1 $bwtParams $strandAlign $toBAM
else 

## Bowtie 2
if [ $optionBowtie == "bwt2" ]; then 

	# Some other options.
	#bwtParams="--time -p3 -k 10 --sensitive-local --no-unal"

	## Samtools parameters
	#Report reads with a mapq of at least 30:
	toBAM="view -q30 -bS"


	################
	# Run function
	echo "running bow1"
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