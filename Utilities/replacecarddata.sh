#! /bin/bash

#
# This script replace the histo_Data of hwwof_0j_shape_8TeV.txt 
# and xwwof_0j_shape_8TeV.txt
# 

export DIR=$1
export NJET=$2
export MH=125

if [ ! $# -eq 2 ]; then
	echo "USAGE: ./replacecarddata.sh njet DIR
	DIR - directory of the card in the smurf card convention $DIR/125/
	NJETS - 0 or 1"
	exit 1
fi


root -l -b -q replacecarddata.C\(\"$DIR/\",$MH,$NJET\)
mv $DIR/$MH/hwwof_${NJET}j_shape_8TeV_temp.txt $DIR/$MH/hwwof_${NJET}j_shape_8TeV.txt 
mv $DIR/$MH/xwwof_${NJET}j_shape_8TeV_temp.txt $DIR/$MH/xwwof_${NJET}j_shape_8TeV.txt 
