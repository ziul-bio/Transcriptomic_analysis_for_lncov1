#!/bin/bash
set -o errexit

# how to run this script 
#bash scripts/run_kallisto_pipeline.sh -experiment PRJNA636804

# Define a list of valid experiment IDs
valid_experiments=("GSE120561" "GSE93028" "PRJNA636804")
experiment=$2
is_valid_experiment=false

# Check if the provided experiment ID is in the valid experiments list
for valid_experiment in "${valid_experiments[@]}"; do
    if [ "$experiment" = "$valid_experiment" ]; then
        is_valid_experiment=true
        break
    fi
done

# If the experiment ID is not valid, exit the script
if [ "$is_valid_experiment" = false ]; then
    echo "Invalid experiment ID: $experiment"
    echo "Valid experiment IDs: ${valid_experiments[*]}"
    exit 1
fi



# Run the pipeline for the specified experiment
if [ "$experiment" = "GSE120561" ]; then
    # Sample - GSE120561 (queen, queenless-worker and worker)
    echo "Running Analysis for GSE120561"
    echo " "
    echo "Doanload the data"
    bash scripts/download_experiments.sh -experiment GSE120561
    echo " "
    echo "Running quantification with kallisto"
    bash scripts/quant.sh -experiment GSE120561
   


elif [ "$experiment" = "GSE93028" ]; then
    # sample experiment id GSE93028 (lncRNA-miRNA-mRNA and reproductive regulation in honey bees)
    echo "Running Analysis for GSE93028"
    echo " "
    echo "Doanload the data"
    bash scripts/download_experiments.sh -experiment GSE93028
    echo " "
    echo "Running quantification with kallisto"
    bash scripts/quant.sh -experiment GSE93028
    


elif [ "$experiment" = "PRJNA636804" ]; then
    # Title: Sex- and caste-specific transcriptomes of larval honey bee (Apis mellifera L.) gonads: DMRT A2 and Hsp83 are differentially expressed and regulated by juvenile hormone
    # Authors: Denyse Cavalcante Lago, Martin Hasselmann, Klaus Hartfelder. 
    echo "Running Analysis for PRJNA636804"
    echo " "
    echo "Doanload the data"
    bash scripts/download_experiments.sh -experiment PRJNA636804
    echo " "
    echo "Running quantification with kallisto"
    bash scripts/quant.sh -experiment PRJNA636804
   
   
fi