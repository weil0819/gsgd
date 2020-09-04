# sudo apt-get update -y
# sudo apt-get install -y gnuplot
set terminal postscript enhanced eps "Helvetica" 38
set size 1,1
set border
set style data histogram
set style histogram clustered gap 1.8
set style fill pattern 1 border -1
set logscale y
set rmargin 1
set ytics
set format y "10^{%L}"
set ylabel "Split Time (s)" offset character -0.2,0 font "Helvetica, 45"
set key left font "Helvetica,32" 
set mytics 5
set ticscale 3
set grid
set xtics right rotate by 45
set xrange [-0.5:3.5]
set xtics font "Helvetica, 35"

set yrange [0.00001:10]
set output '../eps/Eval-II/eval-II-Brightkite.eps'
plot '../dat/eval-II.dat' using ($2)/1000000:xticlabels(1) title "CLUSTER", \
	'' using ($3)/1000000:xticlabel(1) title "MCCNaive", \
	'' using ($4)/1000000:xticlabel(1) title "MCCRandom"

set yrange [0.00001:10]
set output '../eps/Eval-II/eval-II-Gowalla.eps'
plot '../dat/eval-II.dat' using ($5)/1000000:xticlabels(1) title "CLUSTER", \
	'' using ($6)/1000000:xticlabel(1) title "MCCNaive", \
	'' using ($7)/1000000:xticlabel(1) title "MCCRandom"

set yrange [0.00001:10]
set output '../eps/Eval-II/eval-II-Syn1.eps'
plot '../dat/eval-II.dat' using ($8)/1000000:xticlabels(1) title "CLUSTER", \
	'' using ($9)/1000000:xticlabel(1) title "MCCNaive", \
	'' using ($10)/1000000:xticlabel(1) title "MCCRandom"

set yrange [0.00001:100]
set output '../eps/Eval-II/eval-II-Syn2.eps'
plot '../dat/eval-II.dat' using ($11)/1000000:xticlabels(1) title "CLUSTER", \
	'' using ($12)/1000000:xticlabel(1) title "MCCNaive", \
	'' using ($13)/1000000:xticlabel(1) title "MCCRandom"
