#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# HTSeq-count wrapper
# Last modified
# 2017.10.jrm

#################
# For genome mapped sequences

############

## Set parameters for tophat
# cut -f3 TAIR10_withTransposons.gff | sort | uniq # To see the optios for  --type:

display_help() {
	echo "" >&2
	echo "Usage: $0 [option...] " >&2
    echo "Wrapper to run HTSeq on genome mapped RNASeq reads. Needs samtools installed." >&2
    
    echo
    echo "   -h            this useful help"
    echo "   -v            version of bowtie mapped to use [bwt1/bwt2]"
    #echo "   -r            reference [genome/cDNA]"
    #echo "   -a            Path for the annotation file (GFF or GTF file)"
    #echo "---  Optional: "
    #echo "   -p            file with HTSeq arguments
    #		 default: -f sam --stranded yes --mode intersection-strict --idattr ID --type gene"
    #echo
    echo "example: 04a_annotationHTSeq.sh -v bwt1 -a meta/annotation.gff"
    echo ""
    #echo -e "\nA 04_Counts/ directory is created:
   # 	A 'htseqraw' subdirectory that contains the raw output from HTSeq-count.
   # 	A '{bwt1/bwt2}_genome' subdirectory contains the count files.
   # 	In the log folder, a subdirectory 'Counts_{bwt1/bwt2}_genome' with a summary of the annotation." >&2
   # echo ""
    echo -e "\n Assumptions:"
    echo "* A directory '03_alignment/{bwt1,bwt2}_cDNA/' with cDNA/CDS annotated reads in bam format"
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


while getopts ':hv:' option; do
  case "$option" in
    h) display_help
       exit
       ;;
    v) optionBowtie=$OPTARG
    		if [ $optionBowtie == 'bwt1' ] || [ $optionBowtie == 'bwt2' ]
    		    then
    		    echo ""
    		else
    		    echo "Bowtie version not valid please choose [bwt1/bwt2]" >&2
    		    display_help
    		fi
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

# ### If annotation is empty, exit
# if [[ ! -z $annotationPath ]] 
#     then
#     echo ""
# else
# 	echo "Please use a valid reference: -a Path/Reference.{gff,gtf}"
#     display_help
# fi

### If empty, set to default parameters
# if [[  -z  $htseqParams ]] ; then
# 	echo "Bowtie parameters are empty, using default parameters" >&2
	
# 	biotype="gene" #mRNA, CDS, 
# 	attribute="ID" #ID, Parent
# 	mode="intersection-strict" #"union"
# 	htseqParams="-f sam --stranded yes --mode $mode --idattr $attribute --type $biotype"
# else
# 	### Read parameters
# 	htseqParams=`cat $htseqParams`
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

########################################################

## Read
inputFolder=$optionBowtie"_cDNA"

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
# echo `awk '/^__/ { print $0 }' $countFile  > $logCount/$countName.stat.log` #Get all lines starting with double underscore into a new file
#    echo `awk '!/^__/ { print $0 }' $countFile > $outDir/$countName.txt` #Counts file without stats
#done 2>&1 | tee -a $logPath



   
########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: $htseqParams $annotationPath"  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done