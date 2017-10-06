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
	echo "Wrapper to build a reference index using bowtie{1/2}.  " >&2
	echo "If the script detects an existing index it won't create a new one." >&2
	echo "At least one reference (cDNA and/or genome) needed [-c/-g]." >&2
	echo ""
	echo "Usage: $0 [option...] " >&2
    echo ""
    echo
    echo "   -h            this useful help"
    echo "   -v            bowtie version [bwt1/bwt2]"
    echo "   -g            path to the genome fasta"
    echo "   -c            path to the cDNA fasta"
    echo "Optional: "
    echo "   -s            suffix of the reference (default: \"ref\")"
    echo
    echo "example: ./03a_MakeBowtieIndex.sh -v bwt1 -s athIndex -g meta/genome_ath_TAIRv10.fa -c meta/cdna_ath_TAIRv10.fa"
    echo
    #echo "Assumptions:"
    #echo "* A directory with trimmed reads (with extension \"clean\" on a folder named '01_trimmed/')"
    #echo "* The index of the reference is in a 'RefIdx/' directory with format: "
	#echo "  'bwt{v}_reference'. ie, that the directory 'RefIdx/bwt1_cdna/' exists and contains index files"
	#echo "  built with the appropiate bowtie version"
    
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

#In params:
#genomeFasta=$DIR'/meta/genome_ath_TAIRv10.fa'
#cdnaFasta=$DIR'/meta/cdna_ath_TAIRv10.fa'

suffixOut='ref' #Same as BuildIndex script
while getopts ':hv:r:g:c:s:' option; do
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
    g) genomeFasta=$OPTARG
       ;;
    c) cdnaFasta=$OPTARG
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
if [ $optionBowtie == 'bwt1' ] || [ $optionBowtie == 'bwt2' ]
    then
    echo ""
else
	echo "Invalid version please use 'bwt1' or 'bwt2'"
    display_help
fi

if [[  -z  $genomeFasta && -z  $cdnaFasta ]] ; then
	echo "Specify at least one reference fasta file for the index (-g fastaFile or -c fastaFile)"
	display_help
fi



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

## Used to check if index files already exist.
pattern="RefIdx"

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
logBow=$DIR/logs/bowStats"_"$optionAlign$optionBowtie #Create log folder if it doesn't exist
if [ ! -d $logBow ]; then echo `mkdir -p $logBow`; fi
######

## Use the same script's name but add the log extension
## Two ways: split string by dots, keep the first part
# logPath=$DIR/logs/$(basename $BASH_SOURCE | cut -d . -f1).log # What if the script name has more dots?
## Removes specific extension:
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).$optionAlign$optionBowtie.log

## Create a log file and start writing.
echo `touch $logPath` #Create file to fill with log data
echo "Started on `date`" 2>&1 | tee $logPath
############


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

# Out params:
genomeToPath=$DIR"/$pattern/"$optionBowtie"_genome"
cdnaToPath=$DIR"/$pattern/"$optionBowtie"_cdna"


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


   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $bwtParams $strandAlign"  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done