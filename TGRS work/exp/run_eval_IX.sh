#!/bin/bash
# Running time for varying |V|.

for dir in "Syn3" "Syn4" "Syn5" "Syn6" "Syn7" 
do
	for alg in "gdcd" "naive" "random" "greedy"	"swap"
	do
		for alg in "greedy" "swap" "topk"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-IX.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} 0.4 5 50 10  >> ../log/Eval-IX.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-IX.out
		done
	done
done