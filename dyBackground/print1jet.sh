#!/bin/bash

export ANA=$1;

if [ ! $# -eq 1 ]; then
        echo "USAGE: ./print1jet.sh [ANA]
	[ANA] - analysis type, choose from cut, mva and wwxsec"
    exit 1
fi

LUMI=19467 # in pb
export LATEXTABLE=false



echo "| selection | *Rout/in* | N(in)-raw |*N(in) (OF/VZ subt)*| * N(out)-data* | N(out) MC | N(out) data/MC|"

for MASS in 0 115 120 125 130 135 140 145 150 160 170 180 190 200 250 300; do
#for MASS in 0 125; do
    FILE=dyest_mH${MASS}_${LUMI}pb_${ANA}.txt
    if [ $MASS == 0 ]; then
	FILE=dyest_mH${MASS}_${LUMI}pb.txt
    fi
    if [ $ANA == "wwxsec" ]; then
	FILE=dyest_mH${MASS}_${LUMI}pb_wwxsec.txt
    fi
    SELECTION=$MASS
    R=`grep 'R(EE+MM) 1-Jet' $FILE | awk '{print $4""$5""$6""$7""$8}'`
    ZDATA=`grep -C 1 'In-Z peak yield (mm/me/em/ee) after all signal selections 1-Jet' $FILE | tail -n 1 | awk '{print $1"/"$2"/"$3"/"$4}'`
    ZDATASUB=`grep -C 1 'OF/VZ subtracted yields in Z window 1-Jet' $FILE | tail -n 1 | awk '{print $7""$8""$9}'`
    MC=`grep -C 1 'MC estimation in signal region 1-Jet' $FILE | tail -n 1 | awk '{print $7""$8""$9}'`
    DATADRIVEN=`grep -C 1 'data-driven estimate 1-Jet:' $FILE | tail -n 1 | awk '{print $12""$13""$14""$15""$16}'`
    OUTSCALE=`grep -C 1 'data/MC scale factor from Rout/in method 1-Jet' $FILE | tail -n 1 | awk '{print $9""$10""$11}'`
    if [ $LATEXTABLE  == true ]; then 
	echo " $SELECTION \GeV & $ZDATASUB & $R & $DATADRIVEN  & $MC \\\\"
    else
	echo "| $SELECTION | *$R* | $ZDATA | *$ZDATASUB* | *$DATADRIVEN* | $MC | $OUTSCALE |"
    fi
done
