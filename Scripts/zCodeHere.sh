#!/bin/bash
# j rodriguez medina
# Base for
# >> 
# Last modified (yymmdd)
# 


display_help() {
    echo "Usage: $0 [option...] " >&2
    echo
    echo "   -h            this helpful help"
    echo "   -a            here be arguments "
    echo
    # echo some stuff here for the -a or --add-options 
    exit 1
}

# If no argument send help
if [ $# -eq 0 ]
  then
    display_help
fi


while getopts ':ha::' option; do
  case "$option" in
    h) display_help
       exit
       ;;
    a) optionA=$OPTARG
		if (( $optionA > 5 || $optionA < 1 )); then
			echo "Value not valid" >&2
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

### Check that the parameters aren't empty.
if [[  -z  $optionA ]] ; then
	echo "Parameter are empty" >&2
	display_help
fi

################################################################################################
####################################


## Record time
start_time=`date +%s`
## Get absolute path
#Get the full path of the current script: Go back one folder (from scripts to main)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

## Keep track of errors and outputs in a log.
logDir=$DIR/logs #Create log folder if it doesn't exist
if [ ! -d $logDir ]; then echo `mkdir -p $logDir`; fi

# A file for logs
logPath=$DIR/logs/$(basename $BASH_SOURCE .sh).log #Keeps the name of the script, use .log instead of scripts extension
##

echo `touch $logPath` #Create file to fill with log data
##
echo "Started on `date`" 2>&1 | tee $logPath #Start filling 
######



############ Here be functions
#############################################
#Functions help keep code minimal, avoid repeating lines, and just doing the same with diff values. 

#songOftheLonelyFunction () {
  #Write functions like this;
  #Do something with a $variable
  #}

#songOftheLonelyFunction variable1 #Call functions like this
#songOftheLonelyFunction variable2 

############ Code here 
#############################################







########################################################################################################################
########################################################################################################################

## Record time
end_time=`date +%s`

#echo -e "\nParameters used: $bwtParams"  2>&1 | tee -a $logPath
echo -e "\n execution time was `expr $end_time - $start_time` s."  2>&1 | tee -a $logPath
echo -e "\n Done `date`"  2>&1 | tee -a $logPath
##Done