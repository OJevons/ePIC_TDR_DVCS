#! /bin/bash

### Oliver Jevons, University of Glasgow
### Based on files from Stephen Kay, University of York
### 21/09/26
### oliver.jevons@glasgow.ac.uk
### A script to create and submit analysis jobs to the JLab farm.
### As is, takes 3 arguments, campaign, beam energy combination and setting/comment for filelist

RunDir="/w/halla-scshelf2102/sbs/ojevons/eic/ePIC_TDR_DVCS" # Put in the path of your directory here 
echo "Running as ${USER}"
echo "Assuming simulation directory - ${RunDir}"

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
FileList="${RunDir}/bkg_dis/inputFileList_DIS_${Campaign}_${BeamE}_${Setting}.list"
if [ ! -f ${FileList} ]; then
    echo "${FileList} not found!"
    echo "Check path and input arguments carefully!"
    exit 4
fi

echo "Processing - ${FileList}"
echo "Will run as chunks of 10 files (hard coded for now!)"

# Check number of lines in file
NumLines=$(wc -l < "${FileList}")
# Split into 10 line chunks, 3 digit numeric suffix used to identify chunks, 001, 002, etc
split -l 10 -d -a 3 --additional-suffix=.list "${FileList}" "${FileList%.list}_"
# Find number of chunks list was split into and assign to variable to iterate over
NumChunks=$(ls -1 ${FileList%.list}_* 2>/dev/null | wc -l)

Timestamp=$(date +'%d_%m_%Y')
Workflow="ePIC_DISBkg_Analysis_${BeamE}_${Setting}_${USER}_${Timestamp}" # Change this as desired

# Make directories for the input filelists and output ROOT
eval "mkdir ./bkg_dis/${Campaign}_${BeamE}_${Setting}"
eval "mkdir ./rootfiles/${Campaign}_${BeamE}_${Setting}"

# Define a disk space request. Change depending upon your needs. 
Disk_Space=1
for (( i=0; i<$NumChunks; i++ )) # Process all chunks
do
    printf -v ChunkNum "%03d" "${i}" # Convert integer to a padded 3 digit integer
    JOBNAME="DISBkg_Analysis_${Campaign}_${BeamE}_${Setting}_${ChunkNum}"
    # The line below is the "guts" of this script. This is the script or job we will actually run on the farm node. Change the pathing to your script as needed. Submit the args as needed.
    COMMAND="${RunDir}/Process_DIS_Job.sh ${Campaign} ${BeamE} ${Setting}_${ChunkNum}"
    echo "Submitting batch job"
    eval "swif2 add-job -create -workflow ${Workflow} -name ${JOBNAME} -disk 1GB -ram 4GB -account eic -partition production -stdout /farm_out/${USER}/swif/${Workflow}/out/${JOBNAME}.out -stderr /farm_out/${USER}/swif/${Workflow}/err/${JOBNAME}.err ${COMMAND}"
    echo " "
done

eval 'swif2 run ${Workflow}'

exit 0

# Old version making a batch job card and submitting
#batch="${BASE_DIR}/${Campaign}_${BeamE}_${Setting}_FileList_Chunk_${ChunkNum}.txt"
#echo "Running ${batch}"
#cp /dev/null ${batch}
#echo "PROJECT: eic" >> ${batch}
#echo "TRACK: analysis" >> ${batch}    
#batch="${BASE_DIR}/${Campaign}_${BeamE}_${Setting}_FileList_Chunk_${ChunkNum}.txt"
#echo "Running ${batch}"
#cp /dev/null ${batch}
#echo "PROJECT: eic" >> ${batch}
#echo "TRACK: analysis" >> ${batch}    
#echo "MAIL: ${USER}@jlab.org" >> ${batch} # Modify as desired

#eval "swif2 add-jsub ${Workflow} -script ${batch} 2>/dev/null" # Add our created job to the swif2 workflow
#sleep 1
#rm ${batch} # Remove the job script after submission
