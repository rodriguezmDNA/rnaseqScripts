#!/bin/bash
# j rodriguez medina
# Brady lab @ ucdavis
# QC of fastq files (after trimming)
# Last modified
# 2017.07 jrm

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

########################################################################################################################
########################################################################################################################


## Set source/destination folders.
SeqDir="$DIR/RawData/*.fastq.gz"
ToDir=$DIR/00_QC/a_RawQC/
##

echo `mkdir -p $ToDir` #Create directory to store results

###
# Run fastqc
for f in `ls $SeqDir`; do
  echo $f 2>&1 | tee -a $logPath 
  DoQC=`fastqc $f --extract -o $ToDir`   
done >&2 | tee -a $logPath #Only errors. using 2>&1 creates a file too large.


########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

echo -e "\nParameters used: fastqc file --extract -o ToDir" 2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done