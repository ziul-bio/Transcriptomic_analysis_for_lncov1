#!/usr/bin/env bash

# Para terminar a execução do script se algum erro ocorrer
set -o errexit

# how to run:
#bash scripts/quant.sh -experiment PRJNA996350

# Define a list of valid experiment IDs
valid_experiments=("GSE120561" "GSE93028" "PRJNA996350")
experimentID=$2
is_valid_experiment=false


# Check if the provided experiment ID is in the valid experiments list
for valid_experiment in "${valid_experiments[@]}"; do
    if [ "$experimentID" = "$valid_experiment" ]; then
        is_valid_experiment=true
        break
    fi
done


# If the experiment ID is not valid, exit the script
if [ "$is_valid_experiment" = false ]; then
    echo "Invalid experiment ID: $experimentID"
    echo "Valid experiment IDs: ${valid_experiments[*]}"
    exit 1
fi



# Define the paths
KALLISTO_INDEX=data/kallisto_index/index_trans
SEQ_DIR=data/${experimentID}/rawFastq
OUT_DIR=data/${experimentID}/kallisto/quant

mkdir -p ${OUT_DIR}



# Define the file path
FILE="data/kallisto_index/index_trans_AMel"

# Check if the index does not exist
if [ ! -f "$FILE" ]; then
    echo "Index file does NOT exist. Running kallisto index..."
    kallisto index -i ${KALLISTO_INDEX} data/transcripts/AMel_all_transcrits.fasta
else
    echo "Index file found. Running kallisto index..."
fi



# Sample - GSE120561 (queen, queenless-worker and worker)
if [ "$experimentID" = "GSE120561" ]; then

	for SAMPLE in W1 W2 AW1 AW2 Q1 Q2;
	do
		time kallisto quant -t 24 --single --fragment-length=49 --sd=1 -b 1000 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/"$SAMPLE" ${SEQ_DIR}/*_"$SAMPLE"*.fastq.gz
	done
	



# sample experiment id GSE93028 (lncRNA-miRNA-mRNA and reproductive regulation in honey bees)
elif [ "$experimentID" = "GSE93028" ]; then

	for SAMPLE in V01 V02 V03 Q01 Q02 Q03 C01 C02 C03 R01 R02 R03 W01 W02 W03;
	do
		time kallisto quant -t 24 -b 1000 --rf-stranded -i ${KALLISTO_INDEX} -o ${OUT_DIR}/"$SAMPLE" ${SEQ_DIR}/*_"$SAMPLE"_1.fastq.gz ${SEQ_DIR}/*_"$SAMPLE"_2.fastq.gz
	done
	




# sample experiment id PRJNA636804 (Sex- and caste-specific transcriptomes of larval honey bee (Apis mellifera L.) gonads: DMRT A2 and Hsp83 are differentially expressed and regulated by juvenile hormone)
elif [ "$experimentID" = "PRJNA636804" ]; then
	
	for i in Amd1 Amd2 Amd3 Amq1 Amq2 Amq3 Amw1 Amw2 Amw3;
	do
		time kallisto quant -t 24 -b 1000 --rf-stranded -i ${KALLISTO_INDEX} -o ${OUT_DIR}/"$i" ${SEQ_DIR}/$i/*_R1.fastq.gz ${SEQ_DIR}/$i/*_R2.fastq.gz
	done


fi