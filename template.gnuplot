# template.gnuplot
set term png
set output "OUTFILE"
set autoscale
set logscale y
plot "SOURCEFILE" using 1:2 with linespoints