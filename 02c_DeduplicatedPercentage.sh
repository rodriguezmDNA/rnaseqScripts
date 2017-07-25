#!/bin/bash
# j rodriguez medina
# github: @rodriguezmDNA
# Brady lab @ ucdavis
# obtain percentage of deduplicated reads from trimmed libraries. 
# Created
# 2017.07.25_jrm
# Last modified
# 2017.07.25_jrm
###############

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

#### Functions
getData () {
		## Calculate number of reads per library in raw files
	echo "Processing files"  2>&1 | tee -a $logPath
	for each in `ls $trimPath/*$inExt`; do
		Library=`basename $each | cut -d '_' -f 1-6` #Keep the first two fields delimited by an underscore
		#Library=${Library%._fastqc.zip*} #${NAME##*} #Take out basename, aka Path
		Library="$( echo "$Library" | sed -e 's#_fastqc.zip##;' )"
		Total=`unzip -p $each */fastqc_data.txt | grep "Total Deduplicated Percentage"` 
		Total=`echo $Total | cut -d " " -f 4`
		Line=`paste -s -d '\t' <(echo $Library) <(echo $Total)`
		echo $Line;
		echo $Line >> $outTrim
	done >&2 | tee -a $logPath  #Only errors. The output is saved on a file $outRaw 

	}


########################################################################################################################
########################################################################################################################


# Read from
trimPath="$DIR/00_QC/b_trimQC"
inExt=".zip" #Type of file

# Write to 
outPath="$DIR/02_LibrarySizes"
outTrim="$outPath/DeduplicatedPercentage.txt"


# Create folder
echo `mkdir -p $outPath`  2>&1 | tee -a $logPath
# Add header
echo -e "File \t Counts" > $outTrim

getData


#Count also reads that were unassigned to any barcode
#unmatched=`gunzip -dc $DIR/Data/demultiplexed/unmatched.fq.gz | grep "@" | wc -l`
#Line=`paste -s -d '\t' <(echo 'unmatched') <(echo $unmatched)`


# ## Calculate number of reads per library in trimmed files
# echo "Processing trimmed files"  2>&1 | tee -a $logPath
# for each in `ls $trimPath/*$inExt`; do

#         Library=`basename $each | cut -d '_' -f 1-6` #Keep the first two fields delimited by an underscore
#         Library="$( echo "$Library" | sed -e 's#.lane.clean_fastqc.zip##;' )"
#         Total=`unzip -p $each */fastqc_data.txt | grep "Total Sequences"`
#         Total=`echo $Total | cut -d " " -f 3`
#         Line=`paste -s -d '\t' <(echo $Library) <(echo $Total)`
#         echo $Line;
#         echo $Line >> $outTrim

# done 2>&1 | tee -a $logPath #Only redirect errors here 

########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done
