set terminal postscript enhanced eps dashed "NimbusSanL-Regu" 50
set size 2, 0.1
set pointsize 5 
set key below center
set xrange [0:0.1]
set yrange [0:0.1]
unset border
unset tics
set output '../eps/Eval-V/eval-V-title.eps'
plot '../dat/eval-V-Brightkite.dat' using ($1):(($2)/1000) title "Naive" with lp lt 2 lc 7 lw 8 pt 14, \
	'../dat/eval-V-Brightkite.dat' using ($1):(($3)/1000) title "Random" with lp lt 3 lc 7 lw 8 pt 3, \
	'../dat/eval-V-Brightkite.dat' using ($1):(($4)/1000) title "DGCD" with lp lt 4 lc 7 lw 8 pt 6

set terminal postscript enhanced eps dashed "NimbusSanL-Regu" 40
set border
set size 1, 1
set xrange [0.5:4.5]
set xtics ("25" 1, "50" 2, "75" 3, "100" 4)
set logscale y
set ytics
set format y "10^{%L}"
set pointsize 5 
set ylabel "Processing Time (s)" offset character 0,0 font "NimbusRomNo9L-Regu, 46"
set xlabel "{/Symbol g}=" offset character -11,1.6 font "NimbusRomNo9L-Regu, 46"

unset key
set yrange [0.005:100]
set output '../eps/Eval-V/eval-V-Brightkite.eps'
plot '../dat/eval-V-Brightkite.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'../dat/eval-V-Brightkite.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'../dat/eval-V-Brightkite.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

set yrange [0.01:1000]
set output '../eps/Eval-V/eval-V-Gowalla.eps'
plot '../dat/eval-V-Gowalla.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'../dat/eval-V-Gowalla.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'../dat/eval-V-Gowalla.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

set yrange [0.001:100]
set output '../eps/Eval-V/eval-V-Syn1.eps'
plot '../dat/eval-V-Syn1.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'../dat/eval-V-Syn1.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'../dat/eval-V-Syn1.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

set yrange [0.01:1000]
set output '../eps/Eval-V/eval-V-Syn2.eps'
plot '../dat/eval-V-Syn2.dat' using ($1):(($2)/1000000) notitle with lp lt 2 lc 7 lw 8 pt 14, \
	'../dat/eval-V-Syn2.dat' using ($1):(($3)/1000000) notitle with lp lt 3 lc 7 lw 8 pt 3, \
	'../dat/eval-V-Syn2.dat' using ($1):(($4)/1000000) notitle with lp lt 4 lc 7 lw 8 pt 6

