#!/bin/bash

export NJETS=$1;
export INPUTDIR=$2
export OUTPUTDIR=$3
export MVATYPE=$4

if [ ! $# -eq 4 ]; then
    echo "
   USAGE: ./draw2DAll.sh [NJETS] [INPUTDIR] [OUTPUTDIR] [MVATYPE]
    	  
      NJETS - number of jets
      INPUTDIR - location of input historams: (ex) ../cards/2D
      OUTPUTDIR - location to output plots: (ex) plots/2D
      MVATYPE - 2D
    
      (ex) ./draw2DAll.sh 0 ../cards/2D plots/2D 2D"
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
if [ -f ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex ]; then
	rm -f ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
fi

echo writing document ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\documentclass{article}" >  ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\usepackage{times}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\usepackage{epsfig}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "\begin{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex



for MASS in 125 200; do
#for MASS in 110 115 120 125 130 135 140 150 160 180 200 ; do
    
	echo doing $MASS
	#
	# make plots
	#
    #for FLAVOR in of sf; do
    for FLAVOR in of; do
	root -q -l -b histo2d.C+\(\"${INPUTDIR}\/\",${MASS},${NJETS},\"$FLAVOR\",\"${OUTPUTDIR}\/\",\"$MVATYPE\"\);
    done

	#
    # write the figures into a latex version
	#
    
	# data
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#	echo "\epsfig{figure=${OUTPUTDIR}/data_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\epsfig{figure=${OUTPUTDIR}/data_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#	echo "\epsfig{figure=${OUTPUTDIR}/dataerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\epsfig{figure=${OUTPUTDIR}/dataerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\caption{\textbf{Data} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex

	# signal 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/sig_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/sig_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/sigerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/sigerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Signal} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# qqWW 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/qqWW_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/qqWW_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/qqWWerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/qqWWerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{qqWW} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# ggWW 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/ggWW_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/ggWW_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/ggWWerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/ggWWerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{ggWW} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# Wjets 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Wjets_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Wjets_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Wjetserr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Wjetserr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Wjets} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# Top 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Top_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Top_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Toperr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Toperr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Top} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# VV 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/VV_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/VV_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/VVerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/VVerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{VV} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# Zjets 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Zjets_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Zjets_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Zjetserr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Zjetserr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Zjets} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# Wgamma 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Wgamma_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Wgamma_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Wgammaerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Wgammaerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Wgamma} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex 
    
	# Ztt 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Ztt_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Ztt_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/Ztterr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/Ztterr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{Ztt} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex

	# S/B 
	echo "\begin{figure}[htp]" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\begin{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\begin{tabular}{c}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/SoverB_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\epsfig{figure=${OUTPUTDIR}/SoverB_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in} \\\\" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/SoverBerr_${MVATYPE}_mH${MASS}_${NJETS}j_sf.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#    echo "\epsfig{figure=${OUTPUTDIR}/SoverBerr_${MVATYPE}_mH${MASS}_${NJETS}j_of.eps, width=5in}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{tabular}">> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\caption{\textbf{SoverB} : MT vs Mll in SF (left) and OF (right) for mH=${MASS} GeV. Bottom plots are statistical uncertainties relative to the total background yields.}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "\end{center}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "\end{figure}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
    echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex

	echo "\newpage" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
	echo "" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex

done

echo "\end{document}" >> ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex


# latexing all the plots to a pdf file
#echo latexing all the plots to a pdf file at ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.pdf

#latex ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.tex
#dvipdf ${MVATYPE}plots_${NJETS}j.dvi
#mv ${MVATYPE}plots_${NJETS}j.pdf ${OUTPUTDIR}
#rm ${MVATYPE}plots_${NJETS}j.dvi ${MVATYPE}plots_${NJETS}j.log ${MVATYPE}plots_${NJETS}j.aux
#echo plots are located at ${OUTPUTDIR}/${MVATYPE}plots_${NJETS}j.pdf
