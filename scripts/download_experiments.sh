#!/bin/bash

################## Usage ###############
# I used explore SRA to create the downloadable links.
# bash scripts/download_data.sh -experiment PRJNA996350

# Define a list of valid experiment IDs
valid_data=("GSE120561" "GSE93028" "PRJNA996350" PRJNA636804)

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 -experiment <EXPERIMENT_ID>"
    exit 1
fi

# Check if the first argument is "-experiment"
if [ "$1" != "-experiment" ]; then
    echo "Invalid option: $1"
    echo "Usage: $0 -experiment <EXPERIMENT_ID>"
    exit 1
fi

experiment=$2
is_valid_experiment=false

# Check if the provided experiment ID is in the valid data list
for valid_experiment in "${valid_data[@]}"; do
    if [ "$experiment" = "$valid_experiment" ]; then
        is_valid_experiment=true
        break
    fi
done

# If the experiment ID is not valid, exit the script
if [ "$is_valid_experiment" = false ]; then
    echo "Invalid experiment ID: $experiment"
    echo "Valid experiment IDs: ${valid_data[*]}"
    exit 1
fi




########################################################################
if [ "$experiment" = "PRJNA636804" ]; then
########################################################################
    mkdir -p data/PRJNA636804
    # Project ID PRJNA636804
    # Title: Sex- and caste-specific transcriptomes of larval honey bee (Apis mellifera L.) gonads: DMRT A2 and Hsp83 are differentially expressed and regulated by juvenile hormone
    # Authors: Denyse Cavalcante Lago, Martin Hasselmann, Klaus Hartfelder. 

    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/011/SRR11940311/SRR11940311_1.fastq.gz -o data/PRJNA636804/rawFastq/SRR11940311_worker_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/011/SRR11940311/SRR11940311_2.fastq.gz -o data/PRJNA636804/rawFastq/SRR11940311_worker_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/010/SRR11940310/SRR11940310_1.fastq.gz -o data/PRJNA636804/rawFastq/SRR11940310_queen_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/010/SRR11940310/SRR11940310_2.fastq.gz -o data/PRJNA636804/rawFastq/SRR11940310_queen_2.fastq.gz




########################################################################
elif [ "$experiment" = "GSE120561" ]; then
########################################################################
    # Contributor(s)	Duncan EJ, Leask MP, Dearden PK
    # Sample - GSE120561 (queen, queenless-worker and worker)
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120561
    # Sample 4 - GSE120561 (queen, queenless-worker and worker)
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120561
    # Overall design: RNA-seq to measure gene expression in queen,
    # queen right worker and actively laying worker ovaries
    # tissue: ovary
    # Library strategy:	RNA-Seq
    # Library source:	transcriptomic
    # Library selection:	cDNA
    # Instrument model:	Illumina HiSeq 2500
    mkdir -p data/GSE120561/rawFastq

    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/006/SRR7908186/SRR7908186.fastq.gz -o data/GSE120561/SRR7908186_W1_Worker_Pool_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/007/SRR7908187/SRR7908187.fastq.gz -o data/GSE120561/SRR7908187_W2_Worker_Pool_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/008/SRR7908188/SRR7908188.fastq.gz -o data/GSE120561/SRR7908188_AW1_Active_Pool_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/009/SRR7908189/SRR7908189.fastq.gz -o data/GSE120561/SRR7908189_AW2_Active_Pool_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/001/SRR7908191/SRR7908191.fastq.gz -o data/GSE120561/SRR7908191_Q2_Queen_Pool_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/000/SRR7908190/SRR7908190.fastq.gz -o data/GSE120561/SRR7908190_Q1_Queen_Pool_1.fastq.gz

    mkdir -p data/GSE120561
    echo "Downloading GSE120561"
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/006/SRR7908186/SRR7908186.fastq.gz -o data/GSE120561/rawFastq/SRR7908186_W1_Worker_Pool.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/007/SRR7908187/SRR7908187.fastq.gz -o data/GSE120561/rawFastq/SRR7908187_W2_Worker_Pool.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/008/SRR7908188/SRR7908188.fastq.gz -o data/GSE120561/rawFastq/SRR7908188_AW1_Active_Pool.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/009/SRR7908189/SRR7908189.fastq.gz -o data/GSE120561/rawFastq/SRR7908189_AW2_Active_Pool.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/001/SRR7908191/SRR7908191.fastq.gz -o data/GSE120561/rawFastq/SRR7908191_Q2_Queen_Pool.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR790/000/SRR7908190/SRR7908190.fastq.gz -o data/GSE120561/rawFastq/SRR7908190_Q1_Queen_Pool.fastq.gz


########################################################################
elif [ "$experiment" = "GSE93028" ]; then
########################################################################
    # Contributor(s)	Chen X
    # sample experiment id GSE93028 
    # Title:	Integration of lncRNA-miRNA-mRNA reveals novel insights into reproductive regulation in honey bees
    # Organism:	Apis mellifera
    # Experiment type:	Expression profiling by high throughput sequencing Non-coding RNA profiling by high throughput sequencing
    # link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93028

    mkdir -p data/GSE93028/rawFastq
    echo "Downloading GSE93028"

    # (1) ovaries of virgin queens (n=3); Characteristics: tissue: ovaries, age: virgin queen. 
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/008/SRR5136448/SRR5136448_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136448_V01_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/008/SRR5136448/SRR5136448_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136448_V01_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/009/SRR5136449/SRR5136449_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136449_V02_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/009/SRR5136449/SRR5136449_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136449_V02_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/000/SRR5136450/SRR5136450_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136450_V03_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/000/SRR5136450/SRR5136450_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136450_V03_2.fastq.gz

    # (2) ovaries of egg-laying queens (n=3); Characteristics tissue: ovaries, age: normal egg-laying queen
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/001/SRR5136451/SRR5136451_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136451_Q01_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/001/SRR5136451/SRR5136451_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136451_Q01_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/002/SRR5136452/SRR5136452_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136452_Q02_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/002/SRR5136452/SRR5136452_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136452_Q02_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/003/SRR5136453/SRR5136453_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136453_Q03_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/003/SRR5136453/SRR5136453_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136453_Q03_2.fastq.gz

    # (3) ovaries of egg-laying inhibited queens (n=3); Characteristics tissue: ovaries, age: egg-laying inhibited queen treatment: Caged for 7 days
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/004/SRR5136454/SRR5136454_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136454_C01_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/004/SRR5136454/SRR5136454_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136454_C01_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/005/SRR5136455/SRR5136455_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136455_C02_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/005/SRR5136455/SRR5136455_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136455_C02_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/006/SRR5136456/SRR5136456_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136456_C03_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/006/SRR5136456/SRR5136456_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136456_C03_2.fastq.gz

    # (4) ovaries of egg-laying recovered queens (n=3); Characteristics tissue: ovaries, age: egg-laying recovered queen
    # treatment: the queen was released after caged for seven days. And on the eighth day, the queen recovered to normal condition and laid eggs. 
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/007/SRR5136457/SRR5136457_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136457_R01_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/007/SRR5136457/SRR5136457_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136457_R01_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/008/SRR5136458/SRR5136458_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136458_R02_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/008/SRR5136458/SRR5136458_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136458_R02_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/009/SRR5136459/SRR5136459_1.fastq.gz -o data/GSE93028/rawFastq/SRR5136459_R03_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR513/009/SRR5136459/SRR5136459_2.fastq.gz -o data/GSE93028/rawFastq/SRR5136459_R03_2.fastq.gz


    # sample experiment id GSE119256 
    # Title:	Genome-wide analysis of coding and non-coding RNAs in ovary of honey bee workers
    # Organism:	Apis mellifera
    # Source name:	workers_ovary_mRNA including lncRNA
    # Organism:	Apis mellifera
    # Characteristics social caste: honey bee workers
    # tissue: ovaries
    # molecule subtype: mRNA including lncRNA
    # Extracted molecule	total RNA
    # link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119256

    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/009/SRR7770429/SRR7770429_1.fastq.gz -o data/GSE93028/rawFastq/SRR7770429_W01_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/009/SRR7770429/SRR7770429_2.fastq.gz -o data/GSE93028/rawFastq/SRR7770429_W01_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/002/SRR7770432/SRR7770432_1.fastq.gz -o data/GSE93028/rawFastq/SRR7770432_W02_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/002/SRR7770432/SRR7770432_2.fastq.gz -o data/GSE93028/rawFastq/SRR7770432_W02_2.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/005/SRR7770435/SRR7770435_1.fastq.gz -o data/GSE93028/rawFastq/SRR7770435_W03_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR777/005/SRR7770435/SRR7770435_2.fastq.gz -o data/GSE93028/rawFastq/SRR7770435_W03_2.fastq.gz




fi