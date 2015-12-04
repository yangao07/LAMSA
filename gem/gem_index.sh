#!/bin/bash

CMD_PATH=`dirname $0`
GEM_DIR=$CMD_PATH
export PATH=$PATH:$GEM_DIR

if [ $# -ne 1 ]
then
    echo -e usage:
    echo -e $0 ref.fa
    exit 1
fi

#echo "$GEM_DIR/gem-indexer -i $1 -o $1 > /dev/null 2> /dev/null"
$GEM_DIR/gem-indexer -i $1 -o $1 >/dev/null 2>/dev/null
if [[ $? -ne 0 ]];then
    echo -e "\ngem-indexer error. (${1}.log)"
	exit 1
fi
rm ${1}.log
