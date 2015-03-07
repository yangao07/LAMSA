#!/bin/bash
CMD_PATH=`dirname $0`
BWA_DIR=$CMD_PATH/bwa

if [ $# -ne 2 ]
then
    echo -e usage:
    echo -e $0 ref.fa read.fa
    exit 1
fi

#echo -e bwa $1 $2 > ${2}.sai
$BWA_DIR/bwa aln -i 0 $1 $2 > ${2}.sai 2> ${2}.bwa.aln
if [[ $? -ne 0 ]];then
	echo -e bwa-aln error.
	exit 1
fi
#echo -e bwa samse $1 $2.sai $2 > ${2}.bwa.sam
$BWA_DIR/bwa samse $1 ${2}.sai $2 > ${2}.bwa.sam 2>> ${2}.bwa.aln
if [[ $? -ne 0 ]];then
	echo -e bwa-samse error.
	exit 1
fi

rm ${2}.sai
