#!/bin/bash

LUMI=12100 #in pb
LATEXTABLE=false

for ANA in cutbased; do 
    echo "----------------------------------------------------------------------------------------------"
    echo $ANA
    echo "----------------------------------------------------------------------------------------------"
    echo "|  m_H |    0-j meas  |    0-j exp   |     SF      |   1-j meas |    1-j exp   |     SF      |" 
    for MASS in 115 120 130 140 150 160 170 180 190 200; do
	FILE=wwest_mH${MASS}_${LUMI}pb_${ANA}.txt
	DATAZEROJET=`grep 'WW background in signal region in Data (OF) 0-Jet' $FILE | awk '{print $11""$12""$13}'`
	MCZEROJET=`grep 'WW background in signal region in MC (OF) 0-Jet' $FILE | awk '{print $11""$12""$13}'`
	SFZEROJET=`grep 'data/MC scale factor (OF) 0-Jet' $FILE | awk '{print $7""$8""$9}'`
	DATAONEJET=`grep 'WW background in signal region in Data (OF) 1-Jet' $FILE | awk '{print $11""$12""$13}'`
	MCONEJET=`grep 'WW background in signal region in MC (OF) 1-Jet' $FILE | awk '{print $11""$12""$13}'`
	SFONEJET=`grep 'data/MC scale factor (OF) 1-Jet' $FILE | awk '{print $7""$8""$9}'`
	if [ $LATEXTABLE  == true ]; then 
	    echo " $MASS & $SFZEROJET & $SFONEJET \\\\ "
	else
	    echo "|  $MASS |  $DATAZEROJET  |  $MCZEROJET  | $SFZEROJET |  $DATAONEJET  |  $MCONEJET  | $SFONEJET |" 
	fi
    done
    echo "----------------------------------------------------------------------------------------------"
done
