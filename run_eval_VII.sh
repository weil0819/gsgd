#!/bin/bash

for eps in "0.4" "0.6"
do
	for mu in "5" "10"
	do
		for k in "10" "20" "30" "40" "50" 
		do
			for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
			do						
				echo "*********************************************************************************" >> ./out/Eval-VII.out
				# for((i=0;i<3;i++))
				# do
					./gsgd "evalVII" ./data/${dir} ${eps} ${mu} 50 ${k}  >> ./out/Eval-VII.out
				# done
				echo "*********************************************************************************" >> ./out/Eval-VII.out
			done
		done	
	done
done
