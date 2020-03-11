##################### IQ-TREE for nonCoreN3 ###################
currGF=$1
wDirBase="/my/path/to/working/directory"
alnDirBase="$wDirBase/panoct_results"
currCase="nonCoreN3"
wDir="$alnDirBase/$currCase"
phyDir="$wDir/phyDir"
mkdir "$phyDir"
cd "$phyDir"

alnFile="$currGF.fas"

#Run ModelFinder
/my/installation/path/iqtree -s $wDir/$alnFile -m MFP
mv $wDir/$alnFile.* .

#Get model ID
model=$(grep 'Model of substitution' $alnFile.iqtree | perl -pe 's/: /\t/g' | cut -f2)

#Run phylogeny
/my/installation/path/iqtree -s $wDir/$alnFile -m $model -b 100 -nt 5 -redo
mv $wDir/$alnFile.* .

##################################################
