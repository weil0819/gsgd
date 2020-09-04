set terminal postscript enhanced eps dashed "Helvetica" 38
set size 1.2, 1
set border
set style data histogram
set style histogram clustered gap 1.8
set style fill pattern 1 border -1
set logscale y
set rmargin 1
set ytics
set format y "10^{%L}"
set ylabel "Processing Time (s)" offset character -0.5,0 font "Helvetica, 45"
set key left font "Helvetica,32" 
set mytics 5
set ticscale 3
set grid
set format x "%g%%"
set xtics 20
set xrange [-0.5:4.5]
set xtics font "Helvetica, 35"


#unset key
set yrange [0.001:1000000]
set output 'eps/Eval-VIII/eval-VIII-hist.eps'
plot 'data/eval-VIII.dat' using ($2)/1000000:xticlabels(1) title "GDCD", \
	'' using ($3)/1000000:xticlabel(1) title "Naive", \
	'' using ($4)/1000000:xticlabel(1) title "Random", \
	'' using ($5)/1000000:xticlabel(1) title "Greedy", \
	'' using ($6)/1000000:xticlabel(1) title "Swap"