#!/bin/bash

export NJETS=$1;
export INPUTDIR=$2
export OUTPUTDIR=$3
export MVATYPE=$4

if [ ! $# -eq 4 ]; then
        echo "USAGE: ./doAll.sh NJETS INPUTDIR OUTPUTDIR MVATYPE
        NJETS - number of jets 
        INPUTDIR - location of input historams: (ex) ../cards/2D
        OUTPUTDIR - location to output plots
        MVATYPE - mva type, ME or unroll
		
        (ex) ./drawAll.sh 0 ../cards/2D plots/2D/unrolled 2D"
    exit 1
fi


if [ ! -d $INPUTDIR ]; then
    echo Error: Input dir doesnt exist!
    exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
    mkdir -p $OUTPUTDIR
fi


# write the output plots into a tex version
if [ -f ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex ]; then
	rm -f ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
fi

echo writing document ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
echo "\documentclass{article}" >  ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
echo "\usepackage{times}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
echo "\usepackage{epsfig}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
echo "\begin{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex



for MASS in 125; do
    echo doing $MASS
    for FLAVOR in of sf; do
	root -q -l -b mvaoverlay.C+\(\"${INPUTDIR}\/\",${MASS},${NJETS},\"$FLAVOR\",\"${OUTPUTDIR}\/\",\"$MVATYPE\"\);
	#root -q -l -b unrolloverlay.C+\(\"${INPUTDIR}\/\",${MASS},${NJETS},\"$FLAVOR\",\"${OUTPUTDIR}\/\",\"$MVATYPE\"\);

    done

    # write the figures into a latex version
    echo "\begin{figure}[htb]" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\begin{tabular}{cc}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_sf_overlay_lin.eps, width=2.2in}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_of_overlay_lin.eps, width=2.2in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_sf_stack_lin.eps, width=2.2in}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_of_stack_lin.eps, width=2.2in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\caption{${MVATYPE} output (linear version) in SF (left) and OF (right) for mH=${MASS} GeV.}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\begin{tabular}{cc}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_sf_overlay_log.eps, width=2.2in}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_of_overlay_log.eps, width=2.2in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_sf_stack_log.eps, width=2.2in}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/${MVATYPE}_mH${MASS}_${NJETS}j_of_stack_log.eps, width=2.2in}\\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\end{tabular}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\caption{${MVATYPE} output (log version) in SF (left) and OF (right) for mH=${MASS} GeV.}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
done

echo "\end{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex


# latexing all the plots to a pdf file
echo latexing all the plots to a pdf file at ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.pdf

latex ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.tex
dvipdf ${MVATYPE}plots_unroll_${NJETS}j.dvi
mv ${MVATYPE}plots_unroll_${NJETS}j.pdf ${OUTPUTDIR}
rm ${MVATYPE}plots_unroll_${NJETS}j.dvi ${MVATYPE}plots_unroll_${NJETS}j.log ${MVATYPE}plots_unroll_${NJETS}j.aux
echo plots are located at ${OUTPUTDIR}/${MVATYPE}plots_unroll_${NJETS}j.pdf
