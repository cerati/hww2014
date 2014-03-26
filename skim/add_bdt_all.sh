#!/bin/sh

#for MASS in 110 115 120 125 130 135 140 145 150 155 160 170 180 190 200 250 300 350 400 450 500 550 600; do
for MASS in 125; do

	for JETBIN in 0 1; do  

		./add_bdt_53X.sh $JETBIN $MASS /smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/WW/ /smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb_new/WW/ 0  

	done

done

 
