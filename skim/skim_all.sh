#
# script to skim smurf ntuples and add varibles
# ./skim_all.sh /smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/ /smurf/cerati/skims/Run2012_Summer12_SmurfV9_53X-wwxsecfull8tev/
# 

INPUTDIR=$1
OUTPUTDIR=$2

if [ ! $# -eq 2 ]; then
    echo "USAGE: ./skim_all.sh  [INPUTDIR] [OUTPUTDIR]
        INPUTDIR - location of smurf ntuples to skim (e.g. /smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/)
        OUTPUTDIR - location to output skimmed ntuples
		
        (ex) ./skim_all.sh /smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/ /smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb/" 
    exit 1
fi


# check that directories exist

if [ ! -d $INPUTDIR ]; then
        echo Error: Input dir doesnt exist!
        exit 2
fi

if [ ! -d $OUTPUTDIR ]; then
    echo Error: Output dir doesnt exist!
    exit 3
fi

# Do the DY skimming... 
./skim.sh $INPUTDIR $OUTPUTDIR DY 

# Do the WW skimming... 
./skim.sh $INPUTDIR $OUTPUTDIR WW 

# Do the PassFail skimming...
./skim.sh $INPUTDIR $OUTPUTDIR PassFail

# Do the FO skimming...
./skim.sh $INPUTDIR $OUTPUTDIR WGFO 

# Do the SS skimming...
#./skim.sh $INPUTDIR $OUTPUTDIR SS

# do the ZTT control region skimming...
#./skim.sh $INPUTDIR $OUTPUTDIR ZTT

