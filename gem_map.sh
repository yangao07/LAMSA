#!/bin/bash
CMD_PATH=`dirname $0`
GEM_DIR=$CMD_PATH/gem

if [ $# -ne 2 ]
then
	echo -e usage:
	echo -e $0 ref.fa read.fa
	exit 1
fi

ref=/home/tjiang/Mapping_Time_test/Index/GEM/gem_hg19.gem
mismatch=0.04
edit_distance=0.04
indel_length=3

map_out=$2.gem
#$1.gem
$GEM_DIR/gem-mapper -I $ref -i $2 -o $map_out -m $mismatch -e $edit_distance --max-big-indel-length $indel_length --fast-mapping 2> $mapout.log
if [[ $? -ne 0 ]];then
	echo -e gem-mapper error.
	exit 1
fi
