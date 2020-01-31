#!/bin/bash

for eps in "0.4" "0.6"
do
	for mu in "5" "10"
	do
		for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
		do
			for alg in "greedy" "swap" 
			do
				for k in "10" "50"
				do						
					echo "*********************************************************************************" >> ./out/Eval-VII.out
					for((i=0;i<3;i++))
					do
						./gsgd ${alg} ./data/${dir} ${eps} ${mu} 50 ${k}  >> ./out/Eval-VII.out
					done
					echo "*********************************************************************************" >> ./out/Eval-VII.out
				done
			done
		done	
	done
done
