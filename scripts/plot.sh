#!/bin/bash

tmpfile=/tmp/temp.$$
plotDir="$1"
data=${plotDir}/*.plot

if [[ ! -d "plots" ]]; then
	mkdir plots
fi

for plot in $data; do
	dir=${plot%/*}
	graph=${plot%.*}
	graph=${graph##/}
	sed 's,SOURCEFILE,'${plot}',
		 s,OUTFILE,'${graph}'.png,' template.gnuplot > $tmpfile
	gnuplot $tmpfile	
done

rm /tmp/temp.$$