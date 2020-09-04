#!/bin/bash
# Compare clustering and MCC time.

for gamma in "25" "50" "75" "100"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "naive" "random"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-II.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} 0.6 10 ${gamma}  >> ../log/Eval-II.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-II.out
		done
	done
done