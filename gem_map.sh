#!/bin/bash
CMD_PATH=`dirname $0`
GEM_DIR=$CMD_PATH/gem

if [ $# -ne 4 ]
then
	echo -e usage:
	echo -e $0 ref.fa read.fa -d max_matches
	exit 1
fi

ref=/home/tjiang/Mapping_Time_test/Index/GEM/gem_hg19.gem
gem_ref=${1}.gem
mismatch=0.04
edit_distance=0.04
indel_length=3
#max_matches=200 #for l50
max_matches=$4
min_strata=0

map_out=${2}.gem
#$1.gem
#echo "$GEM_DIR/gem-mapper -I $ref -i $2 -o $map_out -m $mismatch -e $edit_distance --max-big-indel-length $indel_length -d $max_matches -D $min_strata --fast-mapping 2> $map_out.log"
$GEM_DIR/gem-mapper -I $ref -i $2 -o $map_out -m $mismatch -e $edit_distance --max-big-indel-length $indel_length -d $max_matches -D $min_strata --fast-mapping 2> $map_out.log
if [[ $? -ne 0 ]];then
	echo -e gem-mapper error.
	exit 1
fi
