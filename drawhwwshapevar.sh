#!/bin/bash

export INPUTDIR=$1
export OUTPUTDIR=$2
export MVATYPE=$3

if [ ! $# -eq 3 ]; then
        echo "USAGE: ./drawhzzshapevar.sh [INPUTDIR] [OUTPUTDIR] [BDT]
        INPUTDIR - location of cards e.g. ../cards/42X_2100pb_MT_shapesyst/
        OUTPUTDIR - location to output plots
        MVATYPE - mvatype, choose from BDT and ME
		
        (ex) ./drawhwwshapevar.sh ../cards/test_19p5ifb/ plots/temp 2D"
     exit 1
fi


if [ ! -d $INPUTDIR ]; then
    echo Error: Input dir doesnt exist!
    exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
    mkdir -p $OUTPUTDIR
fi


#for MH in 110 115 120 125 130 135 140 150 160 180 200 ; do
for MH in 125; do 
    
    root -l -q -b drawhwwshapevar.C+\(\"${MVATYPE}\",${MH},\"${INPUTDIR}\",\"${OUTPUTDIR}\"\);
    
# latexing all the plots to a pdf file
    
    echo latexing all the plots to a pdf file at ${OUTPUTDIR}/shapevar_mH${MH}_${MVATYPE}.tex
    
    latex ${OUTPUTDIR}/shapevar_mH${MH}_${MVATYPE}.tex
    dvipdf shapevar_mH${MH}_${MVATYPE}.dvi
    mv shapevar_mH${MH}_${MVATYPE}.pdf ${OUTPUTDIR}
    rm -f shapevar_mH${MH}_${MVATYPE}.avi shapevar_mH${MH}_${MVATYPE}.aux  shapevar_mH${MH}_${MVATYPE}.log
    echo plots are located at ${OUTPUTDIR}/shapevar_mH${MH}_${MVATYPE}.pdf
    
done
