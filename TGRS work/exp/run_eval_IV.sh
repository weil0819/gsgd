#!/bin/bash
# Running time for varying mu.

for mu in "5" "10" "15" "20"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "naive" "random" "gdcd"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-IV.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} 0.6 ${mu} 50  >> ../log/Eval-IV.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-IV.out
		done
	done
done
