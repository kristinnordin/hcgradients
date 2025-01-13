#!/bin/bash


DESCRIPTION="""
This script is a pipeline for running congrads scripts on functional data in order to calculate connectopic maps.
It supports subject-wise analysis, group-wise analysis, and fitting Trend Surface Models (TSM) on connectopic maps.

Usage:
./run_congrads.sh [options]

Options:
  -h, --help                 Show help message and exit.
  -d, --data_dir <path>      Specify the path to the data directory.
  -s, --subject <list>       Perform subject-wise analysis.
  -g, --group <list>         Perform group-wise analysis.
  -f, --fit_tsm <list>       Fit TSM on connectopic maps.

Examples:
Subject-wise analysis:
./run_congrads.sh -d /path/to/data -s subjects_list.txt
This command will estimate subject-level connectopic maps for the subjects listed in subjects_list.txt.

Group-wise analysis:
./run_congrads.sh -d /path/to/data -g group_list.txt
This command will estimate group-wise connectopic maps for the group of subjects listed in group_list.txt.

Fit TSM on connectopic maps:
./run_congrads.sh -d /path/to/data -f subjects_list.txt
This command will fit TSM models on connectopic maps for the subjects listed in subjects_list.txt.

Expected data structure: 
             
	The 'data directory' is the directory that contains the functional data of each subject. Functional data for each subject should be 
	saved in a separate folder, where the folder name is the subject_ID. The subject_ID is the subject identifier that is listed the the 
	input subjects list. The functional data in the subdirectory is expected to be named as 'func.nii.gz'. Otherwise the script should 
	be modified accordingly.

	ROI mask are expected to be saved in the 'data directory' under a subdirectory named 'masks'. The folder should contain the following
	ROI masks, otherwise the script should be modified accordingly.
		- lh.hippocampus 
		- rh.hippocampus
		- lh.putamen
		- rh.putamen
		- lh.cud_acc
		- rh.cud_acc
		- lh.striatum
		- rh.striatum
		- cortex

Outputs:
	For subject-level analysis, the output will be saved in the same directory as the functional data under a subdirectory named 'congrads_results'.
	For group-level analysis, the output will be saved in a directory named as the group name under the 'data directory'.
	For subject-level connectopic maps estimation and tsm fit, GNU-Parallel is needed for parallel processing of the subjects.

generated: 27 Oct 2021
last modified: 10 Jan 2025 
@farshad.falahati
"""

showHelp() {
    echo "$DESCRIPTION"
    exit 1
}

[ "$1" = "" ] && showHelp

NORM_FLAG=1 # Flag to detemine if normalized connectopic maps should be estimated. Default is to estimate normalized connectopic maps.
SUBJECT_FLAG=0
GROUP_FLAG=0
TSM_FLAG=0

options=$(getopt -l "help,data_dir:,subject:,group:,fit_tsm:" -o "hd:s:g:f:" -a -- "$@")
eval set -- "$options"
while true
do
    case $1 in
        -h|--help)
            showHelp
            exit 0 
            ;;
        -d|--data_dir)
            shift
            export DATA_DIR="$1"
            export MASKS_DIR="$DATA_DIR/masks"
            ;;
        -s|--subject)
            SUBJECT_FLAG=1
            shift
            SUBJECTS_LIST=("$1") 
            ;;
        -g|--group)
            GROUP_FLAG=1
            shift
            GROUP_LIST=("$1") 
            ;;
        -f|--fit_tsm)
            TSM_FLAG=1
            shift
            TSM_LIST=("$1") 
            ;;
        --)
            shift
            break
    esac
    shift
done


# Ensure DATA_DIR is set
if [ -z "$DATA_DIR" ]; then
    echo "Error: data directory is not set. Use -d or --data_dir to specify the data directory."
    exit 1
fi


# Ensure MASKS_DIR exists
if [ ! -d "$MASKS_DIR" ]; then
    echo "Error: 'masks' directory does not exist at $MASKS_DIR. Ensure the 'masks' directory is present within the data directory."
    exit 1
fi


function estimate_subject_maps {
	# Function to estimate connectopic maps for each subject. 
	# Input argument: subject_id that match the subjects in FUNC_DATA
	# directory. 
	
	local START=$(date +%s)
	local SUBJECT=${1}
	local FUNC_DATA="$DATA_DIR/${SUBJECT}/func.nii.gz"
	local RESULTS_DIR="$DATA_DIR/${SUBJECT}/congrads_results"
	mkdir -p $RESULTS_DIR
	
	printf "Loading fMRI data from: $FUNC_DATA \n"
	printf "Saving results to: $RESULTS_DIR \n"

	# List of ROIs
	local ROIS="lh.hippocampus rh.hippocampus lh.putamen rh.putamen lh.cud_acc rh.cud_acc lh.striatum rh.striatum"

	if [ $NORM_FLAG == 1 ]; then 
		# Estimate normalized connectopic maps for the ROIs
		for ROI in $ROIS; do
			./congrads 	-i $FUNC_DATA \
					   	-r $MASKS_DIR/${ROI}.nii.gz \
					  	 -m $MASKS_DIR/cortex.nii.gz \
					   	-o ${RESULTS_DIR}/${ROI}.norm \
					   	-s -z -p -n 3
		done
	else
		# Estimate non-normalized connectopic maps for the ROIs
		for ROI in $ROIS; do
			./congrads 	-i $FUNC_DATA \
					   	-r $MASKS_DIR/${ROI}.nii.gz \
					   	-m $MASKS_DIR/cortex.nii.gz \
					   	-o ${RESULTS_DIR}/${ROI} \
					   	-s -p -n 3
		done
	fi

	local TOTAL_TIME=$(($(date +%s) - $START))
	echo "Subject $SUBJECT completed. Run time: $((TOTAL_TIME/60)) min & $((TOTAL_TIME%60)) sec"
	echo ""

}
export -f estimate_subject_maps


function estimate_group_maps {
	# Function to estimate group-wise connectopic maps (only non-normalized).
	# Input argument: 
	# 	1- list of subjects that connectopic maps should be calculated on them. 
	# Output will be saved in a directory with basename of group  

	local START=$(date +%s)

	printf "\nReading list of subjectsin the group from: \n${1}\n" 

	while read SUBJECT
	do
		local GROUP+="-i ${DATA_DIR}/${SUBJECT}/func.nii.gz "
	done < ${1}
	
	printf "\nGroup input images:\n${GROUP}\n"
	
	local GROUP_NAME=$(basename "${1}" .txt)
	local RESULTS_DIR="${DATA_DIR}/congrads_${GROUP_NAME}"
	mkdir -p $RESULTS_DIR
	printf "\nSaving group-wise analysis results to: \n${RESULTS_DIR}\n"
	
	# List of ROIs
	local ROIS="lh.hippocampus rh.hippocampus lh.putamen rh.putamen lh.cud_acc rh.cud_acc lh.striatum rh.striatum"

	if [ $NORM_FLAG == 1 ]; then 
		# Estimate normalized connectopic maps for the group for each ROI
		for ROI in $ROIS; do
			./congrads 	$GROUP \
					  	-r $MASKS_DIR/${ROI}.nii.gz \
					  	-m $MASKS_DIR/cortex.nii.gz \
					  	-o ${RESULTS_DIR}/${ROI}.norm \
					  	-s -z -p -n 3
		done
	else
		# Estimate non-normalized connectopic maps for the group for each ROI
		for ROI in $ROIS; do
			./congrads 	$GROUP \
					  	-r $MASKS_DIR/${ROI}.nii.gz \
					  	-m $MASKS_DIR/cortex.nii.gz \
					  	-o ${RESULTS_DIR}/${ROI} \
					  	-s -p -n 3
		done
	fi

	local TOTAL_TIME=$(($(date +%s) - $START))
	echo "Group-level analyses completed. Run time: $((TOTAL_TIME/60)) min & $((TOTAL_TIME%60)) sec"
	echo ""
}
export -f estimate_group_maps


function fit_tsm {

	local START=$(date +%s)
	local SUBJECT=${1}
	local RESULTS_DIR="$DATA_DIR/${SUBJECT}"
	mkdir -p $RESULTS_DIR

	# Define the list of ROIs
	ROIS="lh.hippocampus rh.hippocampus lh.putamen rh.putamen lh.cud_acc rh.cud_acc lh.striatum rh.striatum"

	# Fit TSM model of N orders on different cmaps (N = 1, 2, 3, 4, 5)
	for N in {1..5}; do 
		# Loop over ROIs
		for ROI in $ROIS; do
			# select normalised or non-normalised cmaps
			if [ $NORM_FLAG == 1 ]; then 
				local REORDERED_CMAP_FILE="${RESULTS_DIR}/${ROI}.norm/${ROI}.cmaps.rearranged.nii.gz"
				# Check if the cmaps are rearranged then use the rearranged cmaps otherwise use the original cmaps	
				if [ -f "$REORDERED_CMAP_FILE" ]; then
					./congrads 	-i ${RESULTS_DIR}/${ROI}.norm/${ROI}.cmaps.rearranged.nii.gz \
							   	-r ${MASKS_DIR}/${ROI}.nii.gz \
							   	-o ${RESULTS_DIR}/${ROI}.norm/tsm.f${N} \
							   	-F ${N}
				else
					./congrads 	-i ${RESULTS_DIR}/${ROI}.norm/${ROI}.cmaps.nii.gz \
								-r ${MASKS_DIR}/${ROI}.nii.gz \
								-o ${RESULTS_DIR}/${ROI}.norm/tsm.f${N} \
								-F ${N}
				fi
			else
				local REORDERED_CMAP_FILE="${RESULTS_DIR}/${ROI}/${ROI}.cmaps.rearranged.nii.gz"
				# Check if the cmaps are rearranged then use the rearranged cmaps otherwise use the original cmaps	
				if [ -f "$REORDERED_CMAP_FILE" ]; then
					./congrads 	-i ${RESULTS_DIR}/${ROI}/${ROI}.cmaps.rearranged.nii.gz \
							   	-r ${MASKS_DIR}/${ROI}.nii.gz \
							   	-o ${RESULTS_DIR}/${ROI}/tsm.f${N} \
							   	-F ${N}
				else
					./congrads 	-i ${RESULTS_DIR}/${ROI}/${ROI}.cmaps.nii.gz \
								-r ${MASKS_DIR}/${ROI}.nii.gz \
								-o ${RESULTS_DIR}/${ROI}/tsm.f${N} \
								-F ${N}
				fi
			fi
		done
	done 

	local TOTAL_TIME=$(($(date +%s) - $START))
	echo "Subject $SUBJECT completed. Run time: $((TOTAL_TIME/60)) min & $((TOTAL_TIME%60)) sec"
	echo ""
}
export -f fit_tsm


## Run congrads 
printf "=================================================================\n"
printf "Processing started at: $(date)\n"
printf "=================================================================\n"

if [ $SUBJECT_FLAG == 1 ]; then 
	printf "\nEstimation of connectopic maps using congrads\n"
	printf "\nInitiate parallel processing\n"
	mkdir -p $DATA_DIR/parallel_reports
	cat ${SUBJECTS_LIST} | parallel --results $DATA_DIR/parallel_reports --joblog $DATA_DIR/parallel_joblog.txt -j6 --progress estimate_subject_maps
fi

if [ $GROUP_FLAG == 1 ]; then 
	printf "\nEstimation of connectopic maps for a group of subjects\n"
	estimate_group_maps ${GROUP_LIST} 
fi

if [ $TSM_FLAG == 1 ]; then 
	printf "\nFitting Trend Surface Models on Connectopic maps\n"
	printf "\nInitiate parallel processing\n"
	mkdir -p $DATA_DIR/parallel_reports
	cat ${TSM_LIST} | parallel --results $DATA_DIR/parallel_reports --joblog $DATA_DIR/parallel_joblog.txt -j5 --progress fit_tsm
fi 

printf "=================================================================\n"
printf "Processing finished at: $(date)\n"
printf "=================================================================\n"