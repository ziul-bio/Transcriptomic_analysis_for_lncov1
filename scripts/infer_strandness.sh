#!/usr/bin/env bash

# Para terminar a execução do script se algum erro ocorrer
set -o errexit

SAMPLE=GSE120561

KALLISTO_INDEX=data/kallisto_index/index_trans
SEQ_DIR=experiments/${SAMPLE}/rawFastq
OUT_DIR=experiments/${SAMPLE}/kallisto

# Make a subset from a fastqc
zcat ${SEQ_DIR}/SRR7908186_W1_Worker_Pool.fastq.gz | head -n20000 > ${SEQ_DIR}/test.fastq

# Running kallisto quant
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_un ${SEQ_DIR}/test.fastq
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_fr ${SEQ_DIR}/test.fastq --fr-stranded
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_rf ${SEQ_DIR}/test.fastq --rf-stranded

# infering strandness
paste ${OUT_DIR}/test_fr/abundance.tsv ${OUT_DIR}/test_rf/abundance.tsv ${OUT_DIR}/test_un/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > ${OUT_DIR}/test.libtype.txt
cat ${OUT_DIR}/test.libtype.txt | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' >> ${OUT_DIR}/test.libtype.txt




SAMPLE=GSE120561

KALLISTO_INDEX=data/kallisto_index/index_trans_BomTerr
SEQ_DIR=experiments/GSE93274/rawFastq/
OUT_DIR=experiments/GSE93274/kallisto/

# Make a subset from a fastqc
zcat ${SEQ_DIR}/GSM8372213/SRR29668252_OI_rep1.fastq.gz | head -n20000 > ${SEQ_DIR}/test.fastq

# Running kallisto quant
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_un ${SEQ_DIR}/test.fastq
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_fr ${SEQ_DIR}/test.fastq --fr-stranded
kallisto quant --single -l 49 -s 1 -i ${KALLISTO_INDEX} -o ${OUT_DIR}/test_rf ${SEQ_DIR}/test.fastq --rf-stranded

# infering strandness
paste ${OUT_DIR}/test_fr/abundance.tsv ${OUT_DIR}/test_rf/abundance.tsv ${OUT_DIR}/test_un/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > ${OUT_DIR}/test.libtype.txt
cat ${OUT_DIR}/test.libtype.txt | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' >> ${OUT_DIR}/test.libtype.txt
