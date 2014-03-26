#!/bin/bash

# Example:
# ./run_mva.sh 0 130 projlumi datalumi
#
# 1. input files should be available in "data" folder
#  IMPORTANT: this assumes that the data file already inclueds the differential cross-sections for ME method
#    ln -s /smurf/ceballos/tmva/weights weights
# 

# 2. To add the BDT output to the smurfntuples
#    Use 5.28 version root, this is currently available at lxplus

export NJETS=$1;
export MH=$2;
export INPUTDIR=$3
export OUTPUTDIR=$4
export MEFLAG=$5
 
if [ ! $# -eq 5 ]; then
	echo "USAGE: ./add_bdt.sh njets mH inputdir outputdir meflag
        njet - njet bin e.g. 0, 1 or 2
        mH   - SM Higgs mass hypothesis e.g. 120
        inputdir - input directories.. /smurf/yygao/data/EPS/WW/
        outputdir - input directories.. /smurf/yygao/data/EPS/WW/
	meflag - set to 1 if analyzing the inputs from ME code"
	exit 1
fi

# append the inputdir by the number of jets
if [ ${MEFLAG} == "1" ]; then
	INPUTDIR=${INPUTDIR}/${NJETS}j/ME/
else
    INPUTDIR=${INPUTDIR}/${NJETS}j/
fi

### samples must be in "data" folder
rm -f data
ln -s $INPUTDIR data

# check that input directories exist
if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

# make output directory
OUTPUTDIR=${OUTPUTDIR}/mva/$MH/ 
mkdir -p $OUTPUTDIR
rm -f output
ln -s $OUTPUTDIR output

# this is the prefix added to the BDT added ntuples
#export TAG=ntuples2012_PostICHEP_${MH}train_${NJETS}jets; #FIXME : new BDT weights
export TAG=ntuples2012_MultiClass_125train_${NJETS}jets;  
if [ ${NJETS} == "2" ]; then  
    TAG=ntuples2012_PostICHEP_125train_2jets;
fi
export METHODS=BDTG;


# If do the training...
export DO_TRAINING=0;
export trainMVA_smurfFile=trainMVA_smurf.C+;
if [ ${DO_TRAINING} == "1" ]; then
    export SIG_TRAIN=data/hww${MH}.root;
    export BKG_TRAIN=data/qqww${MH}.root;
    if [ ${NJETS} == "2" ]; then
	export SIG_TRAIN=data/hww${MH}.root;
	export BKG_TRAIN=data/top.root;
    fi  
    if [ ! -d $weights ]; then
	mkdir -p weights;
    fi
    ./root-5.28.sh -l -q -b ${trainMVA_smurfFile}\(${NJETS},\"${SIG_TRAIN}\",\"${BKG_TRAIN}\",\"${TAG}\",\"${METHODS}\",${MH}\);
fi


### Fill MVA output information
### MVA output is available in "weights" folder
### an arbitrary list of samples can be added
# define a list of the files to analyze
rm -f list_samples.txt;
cat > list_samples.txt <<EOF
data_ztt
data
qqww
qqww_powheg
ggww
ttbar_powheg
ttbar
tw
wz
zz
www
dyll
wjets
wgamma
zgamma
wgammafo
zgammafo
wglll
wwmcnlo
wwmcnloup
wwmcnlodown
www_PassFail
ttbar_PassFail
data_PassFail
data_ztt_PassFail
qqww_PassFail
qqww_powheg_PassFail
ggww_PassFail
wjets_PassFail
tw_PassFail
wz_PassFail
zz_PassFail
dyll_PassFail
wgamma_PassFail
zgamma_PassFail
wglll_PassFail
ttbar_powheg_PassFail
ttbar_PassFail
wgammafo
zgammafo
data_SS
data_ztt_SS
qqww_SS
qqww_powheg_SS
ggww_SS
wwmcnlo_SS
wwmcnloup_SS
wwmcnlodown_SS
ttbar_powheg_SS
tw_SS
wz_SS
zz_SS
wgamma_SS
zgamma_SS
wgammafo_SS
zgammafo_SS
wglll_SS
dyll_SS
wjets_SS
ttbar_SS
www_SS
hww125_SS
EOF


# ===========================================
# Fill the smurfntuples with the BDT
# ===========================================
export evaluateMVAFile=evaluateMVA_smurf_hww.C+;

# set running perid setting these depend on the lepton efficiency, fakerate and pu histograms 
# check the script directly to make sure all the histograms are read correctly

export PERIOD=3 

## add mva for the HWW signal  
./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/hww${MH}.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",0,$PERIOD\);
mv $OUTPUTDIR/${TAG}_hww${MH}.root $OUTPUTDIR/hww${MH}_${NJETS}j.root

# add mva for the XWW signal at mH=125 GeV 
#if [ $MH -eq 125 ]; then
#    ./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/xww0m125.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",0,$PERIOD\);
#    mv $OUTPUTDIR/${TAG}_xww0m125.root $OUTPUTDIR/xww0m125_${NJETS}j.root
#    ./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/xww0p125.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",0,$PERIOD\);
#    mv $OUTPUTDIR/${TAG}_xww0p125.root $OUTPUTDIR/xww0p125_${NJETS}j.root
#    ./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/xww2p125.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",0,$PERIOD\);
#    mv $OUTPUTDIR/${TAG}_xww2p125.root $OUTPUTDIR/xww2p125_${NJETS}j.root
#fi


# add mva for all background
# add this for 0/1 Jet all mass points and only add 125 for the 2-jet bin
if [ ${NJETS} -lt 2 ] || [ $MH -eq 125 ]; then
    echo njets = $NJETS
    for i in `cat list_samples.txt` ; do
	dataset=${i%%,*};
	echo "filling MVA information in sample: "  $dataset
	PERIOD=3
	echo "period = $PERIOD"
	./root-5.28.sh -l -q -b ${evaluateMVAFile}\(\"data/${dataset}.root\",${MH},\"${METHODS}\",\"${TAG}\",\"\",${NJETS},1,1,\"\",0,$PERIOD\);
	mv $OUTPUTDIR/${TAG}_${dataset}.root $OUTPUTDIR/${dataset}_${NJETS}j.root
    done
fi

rm -f list_samples.txt;

