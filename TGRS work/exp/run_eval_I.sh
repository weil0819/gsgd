#!/bin/bash
# Count number of GSGs.

for eps in "0.4" "0.5" "0.6" "0.7"
do
	for dir in "Brightkite" "Gowalla" "Syn1" "Syn2"
	do
		echo "*********************************************************************************" >> ../log/Eval-I.out
		../../gsgd random ../../data/${dir} ${eps} 5 50  >> ../log/Eval-I.out
		echo "*********************************************************************************" >> ../log/Eval-I.out
	done
done