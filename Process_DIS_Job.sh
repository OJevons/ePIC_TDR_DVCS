#! /bin/bash

### Oliver Jevons, University of Glasgow
### Based on files from Stephen Kay, University of York
### 21/09/26
### oliver.jevons@glasgow.ac.uk
### A script to process DIS file lists with provided arguments.
### As is, takes 3 arguments, campaign, beam energy combination and setting/comment for filelist

RunDir="/w/halla-scshelf2102/sbs/ojevons/eic/ePIC_TDR_DVCS" # Put in the path of your directory here 
# Put in the path of your directory here (where your eic-shell is)
echo "Running as ${USER}"
echo "Assuming directory - ${RunDir}"

Campaign=$1 # First arg is the campaign version to run
if [[ -z "$1" ]]; then
    echo "I need a campaign version to run!"
    echo "Please provide a campaign version to run as the first argument, this should map to an available file list."
    echo "E.g. 26_04_1"
    exit 2
fi
BeamE=$2 # Second argument is a beam energy combination
if [[ -z "$2" ]]; then
    echo "I need a beam energy combination to analyse!"
    echo "Please provide a beam energy combination as the second argument, this should map to an available file list."
    echo "E.g. 10x130"
    exit 3
fi
Setting="${3:-}" # Assigns third argument if it was provided, set to blank if not

# Construct file list from arguments and check it exists
FileList="${RunDir}/bkg_dis/${Campaign}_${BeamE}_${Setting}/inputFileList_DIS_${Campaign}_${BeamE}_${Setting}.list"
if [ ! -f ${FileList} ]; then
    echo "${FileList} not found!"
    echo "Check path and input arguments carefully!"
    exit 4
fi

echo "Processing - ${FileList}"

export EICSHELL=/w/halla-scshelf2102/sbs/${USER}/eic/eic-shell # Must point to where your eic-shell is!
cat <<EOF | $EICSHELL
cd ${RunDir}
source setup.sh
root -l -b -q 'run_ePIC_DISBkg.C("${Campaign}", "${BeamE}", "${Setting}")'
EOF

sleep 2

#rm $FileList

exit 0
