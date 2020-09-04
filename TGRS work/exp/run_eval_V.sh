#!/bin/bash
# Running time for varying gamma.

for gamma in "25" "50" "75" "100"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		for alg in "naive" "random" "gdcd"
		do		 
			echo "*********************************************************************************" >> ../log/Eval-V.out
			for((i=0;i<3;i++))
			do
				../../gsgd ${alg} ../../data/${dir} 0.6 10 ${gamma}  >> ../log/Eval-V.out
			done			
			echo "*********************************************************************************" >> ../log/Eval-V.out
		done
	done
done