#!/bin/bash


benchmark="$1"
fn="../benchmark_results/${benchmark}_`hg id -b`_`hg id -i`"
filename="${fn}.txt"

fnr=2
while [ -f $filename ]
do
	filename="${fn}_${fnr}.txt"
	((fnr++))
done

[ -d .buildO ] && rm -r .buildO
[ -f NetworKit-Tests-O ] && NetworKit-Tests-O
scons --optimize=O --target=Tests -j 64

./NetworKit-Tests-O --tests --loglevel=ERROR --gtest_filter=*${benchmark}* > $filename
