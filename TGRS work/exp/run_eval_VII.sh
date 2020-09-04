#!/bin/bash
# Running time for varying k.

for k in "10" "20" "30" "40" "50" 
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "greedy" "swap" "topk"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-VII.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} 0.6 5 50 ${k}  >> ../log/Eval-VII.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-VII.out
		done
	done
done