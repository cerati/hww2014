#!/bin/bash

# publish ww preselection level plots 

export INPUTDIR=$1
export OUTPUTDIR=$2

if [ ! $# -eq 2 ]; then
        echo "USAGE: ./publish.sh INPUTDIR OUTPUTDIR
        INPUTDIR  - location of the plots to be published, 
        OUTPUTDIR - location to publish the plots"
    exit 1
fi

for NJETS in 0 1 2; do 
    for FLAVOR in ee mm em me sf of; do 
	export output=${OUTPUTDIR}/${NJETS}j/${FLAVOR}/
	mkdir -p $output
	cp $INPUTDIR/hww_analysis16_0_ALL_${FLAVOR}_${NJETS}j_*.png  $output
	if [ ! -d $output ]; then
	    mkdir -p $output
	fi
	~/publishplot.pl ${OUTPUTDIR}/${NJETS}j/${FLAVOR}/ --t "${FLAVOR} ${NJETS}-Jet in ww level"
    done
done