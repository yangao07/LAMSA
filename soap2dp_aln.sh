#!/bin/bash

SOAP2DP_DIR=./soap2-dp
CUR_DIR=$PWD

if [ $# -ne 2 ]
then
    echo -e usage:
    echo -e $0 ref.fa read.fa
    exit 1
fi

# get relative path for soap2dp
cd $SOAP2DP_DIR
ABS_SOAP=$PWD

N_L=`echo $ABS_SOAP | tr -cd [^/] | wc -c`
echo $N_L

cd $CUR_DIR
ABS_REF_DIR=`dirname $1`
cd $ABS_REF_DIR
ABS_REF=$PWD/`basename $1`

REL_REF=..$ABS_REF
for (( i=1;i<$N_L;++i ))
do
    REL_REF="../"$REL_REF
done

cd $CUR_DIR
ABS_READ_DIR=`dirname $2`
cd $ABS_READ_DIR
ABS_READ=$PWD/`basename $2`


cd $ABS_SOAP 
echo "./soap2-dp single ${REL_REF}.index ${ABS_READ} -h 2 -m 3e -b 2 > ${ABS_READ}.soap2dp.aln"
./soap2-dp single -h 2 -m 3e -b 2 ${REL_REF}.index ${ABS_READ} > ${ABS_READ}.soap2dp.aln


cd $CUR_DIR

