#!/bin/bash
BWA_DIR=./bwa

if [ $# -ne 1 ]
then
    echo -e usage:
    echo -e $0 ref.fa
    exit 1
fi

#echo -e bwa $1 $2 > ${2}.sai
$BWA_DIR/bwa index $1 2> ${1}.bwa.index.build
if [[ $? -ne 0 ]];then
	echo -e bwa-index error.
	exit 1
fi
