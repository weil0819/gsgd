#!/bin/bash
# Running time for varying epsilon.

for eps in "0.4" "0.5" "0.6" "0.7"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "naive" "random" "gdcd"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-III.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} ${eps} 10 50  >> ../log/Eval-III.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-III.out
		done
	done
done
