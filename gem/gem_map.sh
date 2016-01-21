#!/bin/bash
# argv: 
#      $1: ref_prefix
#      $2: read_prefix
#      $3: mis_rate
#      $4: edit_distance
#      $5: min_match_bases
#      $6: max_match_per_read
#      $7: thread_n

CMD_PATH=`dirname $0`
GEM_DIR=$CMD_PATH

#if [ $# -ne 7 ]
#then
#	echo -e usage:
#	echo -e $0 ref.fa read.fa mis_rate edit_distance min_match_bases max_matches_per_read thread_n
#	exit 1
#fi

ref=${1}.gem
read_f=$2
mismatch=$3
edit_distance=$4
indel_length=3
min_match_bases=$5
#max_matches=200
max_matches=$6
min_strata=0
thread_n=$7

map_out=${2}.gem
#echo "$GEM_DIR/gem-mapper -I $ref -i $read_f -o $map_out -m $mismatch -e $edit_distance --min-matched-bases $min_match_bases --max-big-indel-length $indel_length -d $max_matches -D $min_strata -T $thread_n --fast-mapping 2> $map_out.log"
time $GEM_DIR/gem-mapper -I $ref -i $read_f -o $map_out -m $mismatch -e $edit_distance --min-matched-bases $min_match_bases --max-big-indel-length $indel_length -d $max_matches -D $min_strata -T $thread_n --fast-mapping 2> $map_out.log
if [[ $? -ne 0 ]];then
	exit 1
fi
