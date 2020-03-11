# Create outDirs & outfiles for likelihood ratio tests (LRTs) 
# Calculate LTRs with the appropriate R script
# Get lists of gene families under positive selection and count results
# Map gene family to annotation table (one of panoct's input files)

# For core gene families
wDirBase="/my/path/to/working/directory"
alnDirBase="$wDirBase/panoct_results"
currCase="coreProts"
wDir="$alnDirBase/$currCase"
outDir="$wDir/codeml_selection"
mkdir "$outDir"
cd "$outDir"

echo -e "contrast\tpval\tH0\tH1\trejectH0\tgeneFam" > LRTs_001.tab
echo -e "contrast\tpval\tH0\tH1\trejectH0\tgeneFam" > LRTs_05.tab
Rscript calc_LRTs.R $currCase $alnDirBase/ $outDir/

cat LRTs_001.tab | grep yes | cut -f6 > positiveSelFams_coreProts_001.list
cat LRTs_05.tab | grep yes | cut -f6 > positiveSelFams_coreProts_05.list
wc -l positiveSelFams_coreProts_001.list
wc -l positiveSelFams_coreProts_05.list 

inFile="$wDir/noREC.$currCase.gf.list"
annotFile="$alnDirBase/geneFamilies_panoct.txt"
awk -F '\t' 'FNR==NR {hash[$0]; next} $1 in hash' $inFile $annotFile > geneFamilies_panoct.coreProts.txt


# For non-core gene families
wDirBase="/my/path/to/working/directory"
alnDirBase="$wDirBase/panoct_results"
currCase="nonCoreN3"
wDir="$alnDirBase/$currCase"
outDir="$wDir/codeml_selection"
mkdir "$outDir"
cd "$outDir"
echo -e "contrast\tpval\tH0\tH1\trejectH0\tgeneFam" > LRTs_001.tab
echo -e "contrast\tpval\tH0\tH1\trejectH0\tgeneFam" > LRTs_05.tab
Rscript calc_LRTs.R $currCase $alnDirBase/ $outDir/
cat LRTs_001.tab | grep yes | cut -f6 > positiveSelFams_nonCoreN3_001.list
cat LRTs_05.tab | grep yes | cut -f6 > positiveSelFams_nonCoreN3_05.list
wc -l positiveSelFams_nonCoreN3_001.list
wc -l positiveSelFams_nonCoreN3_05.list

inFile="$wDir/noREC.$currCase.gf.list"
annotFile="$alnDirBase/geneFamilies_panoct.txt"
awk -F '\t' 'FNR==NR {hash[$0]; next} $1 in hash' $inFile $annotFile > geneFamilies_panoct.nonCoreN3.txt


