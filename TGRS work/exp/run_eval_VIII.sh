#!/bin/bash
# Running time for varying epsilon.

for eps in "0.4" "0.5" "0.6" "0.7"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "greedy" "swap" "topk"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-VIII.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} ${eps} 5 50 10 >> ../log/Eval-VIII.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-VIII.out
		done
	done
done