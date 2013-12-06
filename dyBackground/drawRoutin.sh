#!/bin/bash


LUMI=12100

for NJETS in 0; do

# write the output plots into a tex version
    if [ -f R${NJETS}j_${LUMI}pb.tex ]; then
        rm -f R${NJETS}j_${LUMI}pb.tex
    fi
    
    echo writing document R${NJETS}j.tex
    
    echo "\documentclass{article}" > R${NJETS}j_${LUMI}pb.tex
    echo "\usepackage{times}" >> R${NJETS}j_${LUMI}pb.tex
    echo "\usepackage{epsfig}" >> R${NJETS}j_${LUMI}pb.tex
    echo "\begin{document}" >> R${NJETS}j_${LUMI}pb.tex
    
    #for MASS in 0 120 125 130 135 140 150 160 170 180 190 200; do
	for MASS in 0; do 
	echo "\begin{figure}[htb]" >> R${NJETS}j_${LUMI}pb.tex
	echo "\begin{center}" >> R${NJETS}j_${LUMI}pb.tex
	echo "\begin{tabular}{ccc}" >> R${NJETS}j_${LUMI}pb.tex
	
	for SAMPLE in dy; do  
	    echo doing $MASS for sample $SAMPLE
	    echo doing "root -l -q draw.C\($MASS,$NJETS,$LUMI,\"$SAMPLE\"\);"
	    root -b -l -q draw.C++\($MASS,$NJETS,$LUMI,\"$SAMPLE\"\);
	    echo "\epsfig{figure=Routin_ee_${NJETS}Jet_mH${MASS}_${LUMI}pb_${SAMPLE}.eps, width=1.6in}" >> R${NJETS}j_${LUMI}pb.tex
	    echo "\epsfig{figure=Routin_mm_${NJETS}Jet_mH${MASS}_${LUMI}pb_${SAMPLE}.eps, width=1.6in}" >> R${NJETS}j_${LUMI}pb.tex
	    echo "\epsfig{figure=Routin_${NJETS}Jet_mH${MASS}_${LUMI}pb_${SAMPLE}.eps, width=1.6in}" >> R${NJETS}j_${LUMI}pb.tex
	done
	
	echo "\end{tabular}" >> R${NJETS}j_${LUMI}pb.tex
	echo "\caption{ Routin dependence on the met for ee (left) mm (middle) and ee/mm for mH=${MASS} GeV.}" >> R${NJETS}j_${LUMI}pb.tex
	echo "\end{center}" >> R${NJETS}j_${LUMI}pb.tex
	echo "\end{figure}" >> R${NJETS}j_${LUMI}pb.tex
	
	
    done
    
    
    echo "\end{document}" >> R${NJETS}j_${LUMI}pb.tex
    
    latex R${NJETS}j_${LUMI}pb.tex
    dvipdf R${NJETS}j_${LUMI}pb.dvi
# clean up
    rm -f R${NJETS}j_${LUMI}pb.dvi
    rm -f R${NJETS}j_${LUMI}pb.log
    rm -f R${NJETS}j_${LUMI}pb.aux
    
done
