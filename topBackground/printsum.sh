#!/bin/bash

# this scripts outputs the table for the mass-dependent VBF top estimations

LUMI=3553

echo ""
for MASS in 110 115 120 125 130 135 140 145 150 155 160 170 180 190 200 250 300 350 400 450 500 550 600; do
    FILE=topest_${LUMI}pb_mH${MASS}.txt
    MC=`grep 'Estimated top events in simulation (SF+OF) 2-Jet' $FILE | awk '{print $8""$9""$10}'`
    DATA=`grep 'Estimated top events in data (SF+OF) 2-Jet' $FILE | awk '{print $8""$9""$10}'`
    echo " $MASS \GeV & $MC & $DATA  \\\\ "
done
echo "----------------------------------------------------------------------------------------------"
