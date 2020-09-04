set terminal postscript enhanced eps dashed "NimbusSanL-Regu" 50
set size 2, 0.1
set pointsize 5 
set key below center
set xrange [0:0.1]
set yrange [0:0.1]
unset border
unset tics
set output 'eps/Eval-VII/evalVII-title.eps'
plot 'data/k-Brightkite.dat' using ($1):(($2)/1000) title "TopK" with lp lt 2 lc 7 lw 8 pt 14, \
	'data/k-Brightkite.dat' using ($1):(($3)/1000) title "Greedy" with lp lt 3 lc 7 lw 8 pt 3, \
	'data/k-Brightkite.dat' using ($1):(($4)/1000) title "Swap" with lp lt 4 lc 7 lw 8 pt 6

set terminal postscript enhanced eps dashed "NimbusSanL-Regu" 40
set border
set size 1, 1
set xrange [0.5:5.5]
set xtics ("10" 1, "20" 2, "30" 3, "40" 4, "50" 5)
set logscale y
set ytics
set format y "10^{%L}"
set pointsize 5 
set ylabel "Processing Time (s)" offset character 0,0 font "NimbusRomNo9L-Regu, 46"
set xlabel "k=" offset character -11,1.6 font "NimbusRomNo9L-Regu, 46"

unset key
set yrange [0.01:1]
set output 'eps/Eval-VII/evalVII-k-Brightkite.eps'
plot 'data/k-Brightkite.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'data/k-Brightkite.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'data/k-Brightkite.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

set yrange [0.01:10]
set output 'eps/Eval-VII/evalVII-k-Gowalla.eps'
plot 'data/k-Gowalla.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'data/k-Gowalla.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'data/k-Gowalla.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

set yrange [0.001:0.1]
set output 'eps/Eval-VII/evalVII-k-Syn1.eps'
plot 'data/k-Syn1.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'data/k-Syn1.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'data/k-Syn1.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6	

set yrange [0.01:10]
set output 'eps/Eval-VII/evalVII-k-Syn2.eps'
plot 'data/k-Syn2.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'data/k-Syn2.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'data/k-Syn2.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

