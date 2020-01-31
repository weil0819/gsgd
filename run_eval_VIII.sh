#!/bin/bash

for eps in "0.6"
do
	for mu in "10"
	do
		for dir in "Syn3" "Syn4" "Syn5" "Syn6" "Syn7"
		do	
			for alg in "gdcd" "naive" "random" "greedy"	"swap"		
			do	
				echo "*********************************************************************************" >> ./out/Eval-VIII.out
					./gsgd ${alg} ./data/${dir} ${eps} ${mu} 50 10  >> ./out/Eval-VIII.out
				echo "*********************************************************************************" >> ./out/Eval-VIII.out
			done
		done	
	done
done