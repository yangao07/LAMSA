#!/bin/bash
# argv: 
#      $1: ref_prefix
#      $2: read_prefix
#      $3: max_match_per_read
#      $4: thread_n

CMD_PATH=`dirname $0`
GEM_DIR=$CMD_PATH

if [ $# -ne 4 ]
then
	echo -e usage:
	echo -e $0 ref.fa read.fa -d max_matches
	exit 1
fi

ref=${1}.gem
read_f=$2
mismatch=0.04
edit_distance=0.04
indel_length=3
#max_matches=200
max_matches=$3
min_strata=0
thread_n=$4

map_out=${2}.gem
#echo "$GEM_DIR/gem-mapper -I $ref -i $read_f -o $map_out -m $mismatch -e $edit_distance --max-big-indel-length $indel_length -d $max_matches -D $min_strata -T $thread_n --fast-mapping 2> $map_out.log"
$GEM_DIR/gem-mapper -I $ref -i $read_f -o $map_out -m $mismatch -e $edit_distance --max-big-indel-length $indel_length -d $max_matches -D $min_strata -T $thread_n --fast-mapping 2> $map_out.log
if [[ $? -ne 0 ]];then
	echo -e gem-mapper error.
	exit 1
fi
