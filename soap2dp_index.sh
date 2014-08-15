#!/bin/bash
SOAP2DP_DIR=./soap2-dp

if [ $# -ne 1 ]
then
    echo -e usage:
    echo -e $0 ref.fa
    exit 1
fi

$SOAP2DP_DIR/soap2-dp-builder $1 2> ${1}.soap2dp.build
if [[ $? -ne 0 ]];then
	echo -e soap2-dp-build error.
	exit 1
fi
