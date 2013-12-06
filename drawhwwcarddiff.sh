#!/bin/bash

export NJETS=$1;
export MVATYPE=$2
export INPUTDIR1=$3
export INPUTDIR2=$4
export INPUTNAME1=$5
export INPUTNAME2=$6
export OUTPUTDIR=$7

if [ ! $# -eq 7 ]; then
        echo "USAGE: ./drawhwwcarddiff.sh NJETS MVATYPE INPUTDIR1 INPUTDIR2 INPUTNAME1 INPUTNAME2 OUTPUTDIR
        NJETS - number of jets 
        MVATYPE - mva type, ME, BDT or 2D
        INPUTDIR1 - location of first set of cards
        INPUTDIR2 - location of second set of cards
        INPUTNAME1 - location of first set of cards
        INPUTNAME2 - location of second set of cards
        OUTPUTDIR - location to output plots"
    exit 1
fi


if [ ! -d $INPUTDIR1 ]; then
    echo Error: Input dir 1 doesnt exist!
    exit 2
fi

if [ ! -d $INPUTDIR2 ]; then
    echo Error: Input dir 2 doesnt exist!
    exit 2
fi


if [ ! -d $OUTPUTDIR ]; then
    mkdir -p $OUTPUTDIR
fi


# write the output plots into a tex version
if [ -f ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex ]; then
	rm -f ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
fi

echo writing document ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\documentclass{article}" >  ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\usepackage{times}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\usepackage{epsfig}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\begin{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex



#for MASS in 110 115 120 125 130 135 140 145 150 160 170 180 190 200 250 300 350 400 450 500 550 600; do
for MASS in 125; do
    echo doing $MASS
    for FLAVOR in of; do
	root -q -l -b drawhwwcarddiff.C+\(\"${INPUTDIR1}\/\",\"${INPUTDIR2}\/\",${MASS},${NJETS},\"$FLAVOR\",\"${OUTPUTDIR}\/\",\"$MVATYPE\",\"$INPUTNAME1\",\"$INPUTNAME2\"\);
    # write the figures into a latex version
    echo "\begin{figure}[htb]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{tabular}{cc}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_ggH_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_qqH_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_qqWW_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_ggWW_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Top_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Wjets_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_VV_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Wgamma_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_qqWW_CMS_hww_MVAWWBoundingUp_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_qqWW_CMS_hww_MVAWWBoundingDown_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Top_CMS_hww_MVATopBoundingUp_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    if [ $FLAVOR == "sf" ]; then
      echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Zjets_CMS_hwwsf_${NJETS}j_MVAZBoundingUp_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    fi
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Data_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Top_CMS_hww_MVATopBoundingDown_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    
    #echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_histo_Ztt_mH${MASS}_${NJETS}j_${FLAVOR}.eps, width=2.0in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{${MVATYPE} output differences for \bf{${FLAVOR} ${NJETS} jet bin} of mH=${MASS} GeV.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\clearpage" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    done
done

echo "\end{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex


# latexing all the plots to a pdf file
echo latexing all the plots to a pdf file at ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.pdf

latex ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
dvipdf ${MVATYPE}plots_${NJETS}j.dvi
mv ${MVATYPE}plots_${NJETS}j.pdf ${OUTPUTDIR}
echo plots are located at ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.pdf
