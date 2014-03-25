
#
# script to skim smurf ntuples
# choose from WW preselection, and PassFail Selections
# 

INPUTDIR=$1
OUTPUTDIR=$2
SELECTION=$3

if [ ! $# -eq 3 ]; then
    echo "USAGE: ./skim.sh [INPUTDIR] [OUTPUTDIR] [skim]
        INPUTDIR - location of smurf ntuples to skim (e.g./smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/ )
        OUTPUTDIR - location to output skimmed ntuples
        SELECTION - selection, choose from WW, PassFail, PassFailSS, ZTT
        
        (ex) ./skim_all.sh /smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/ /smurf/jaehyeok/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets_19p5ifb/ WW" 

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

# loop over root files in input dir
# and do the skim root script

if [ "$SELECTION" == 'WW' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data_ztt.root
data.root
data_3l.root
ttbar.root
tw.root
qqww.root
qqww_powheg.root
ggww.root
zz.root
wz.root
wjets.root
dyll.root
wgamma.root
zgamma.root
wglll.root
www.root
wwmcnlo.root
wwmcnloup.root
wwmcnlodown.root
ttbar_powheg.root
hww110.root
hww115.root
hww120.root
hww125.root
hww130.root
hww135.root
hww140.root
hww145.root
hww150.root
hww155.root
hww160.root
hww170.root
hww180.root
hww190.root
hww200.root
hww250.root
hww300.root
hww350.root
hww400.root
hww450.root
hww500.root
hww550.root
hww600.root 
xww0m125.root  
xww0p125.root  
xww2p125.root 
vhtt110.root  
vhtt115.root  
vhtt120.root  
vhtt125.root  
vhtt130.root  
vhtt135.root  
vhtt140.root  
vhtt145.root  
vhtt150.root  
vhtt155.root  
vhtt160.root  
EOF
fi



if [ "$SELECTION" == 'PassFail' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data.root
qqww.root
qqww_powheg.root
ggww.root
ttbar_powheg.root
tw.root
wz.root
zz.root
wgamma.root
zgamma.root
wglll.root
dyll.root
wjets.root
ttbar.root
www.root
wgammafo.root
zgammafo.root
data_ztt.root
EOF
fi


if [ "$SELECTION" == 'SS' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data_ztt.root
data.root
qqww.root
qqww_powheg.root
ggww.root
ttbar_powheg.root
tw.root
wz.root
zz.root
wgamma.root
zgamma.root
wglll.root
dyll.root
wjets.root
ttbar.root
www.root
wgammafo.root
zgammafo.root
wwmcnlo.root 
wwmcnloup.root
wwmcnlodown.root
EOF
fi

if [ "$SELECTION" == 'WGFO' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
wgammafo.root
zgammafo.root
EOF
fi

if [ "$SELECTION" == 'ZTT' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data_ztt.root
data.root
qqww.root
qqww_powheg.root
ggww.root
ttbar_powheg.root
tw.root
wz.root
zz.root
wgamma.root
zgamma.root
wglll.root
dyll.root
wjets.root
ttbar.root
www.root
wgammafo.root
zgammafo.root
wwmcnlo.root 
wwmcnloup.root
wwmcnlodown.root
EOF
fi

if [ "$SELECTION" == 'DY' ]; then
rm -f list_samples.txt
cat > list_samples.txt <<EOF
data.root
dyll.root
wz.root
zz.root
EOF
fi

# Do the skimming...
for FILE in `cat list_samples.txt` ; do

    for JETBIN in 0 1 2 ; do 
	outputdir=$OUTPUTDIR/$SELECTION/${JETBIN}j/
	if [ "$SELECTION" == 'PassFail' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
	if [ "$SELECTION" == 'WGFO' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
	if [ "$SELECTION" == 'SS' ]; then
	    outputdir=$OUTPUTDIR/WW/${JETBIN}j/
	fi
    if [ "$SELECTION" == 'ZTT' ]; then
        outputdir=$OUTPUTDIR/WW/${JETBIN}j/
    fi

    if [ "$SELECTION" == 'DY' ]; then 
        if [ $JETBIN == 0 ]; then 
        mkdir -p $OUTPUTDIR/dyskim/
	    echo doing "root -l -b -q dyskim.C+\(\"$INPUTDIR\",\"$FILE\",\"$OUTPUTDIR/dyskim/\"\);"
	    root -l -b -q dyskim.C+\(\"$INPUTDIR\",\"$FILE\",\"$OUTPUTDIR/dyskim/\"\);
        fi
    else 
	    mkdir -p $outputdir
	    echo doing "root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);"
	    root -l -b -q smurfproducer.C+\(\"$INPUTDIR\",\"$FILE\",\"$outputdir\",\"$SELECTION\",$JETBIN\);
    fi
    
    done

done  


