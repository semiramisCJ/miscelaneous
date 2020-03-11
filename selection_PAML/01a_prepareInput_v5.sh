# *Note: panoctDir has the results of PanOCT's run; for an example of how to obtain them,
# check the following directory in the same repository
# ../pangenomes-core_phylogenetic_trees/02_run_panoct.sh

# This script shows how to get the corresponding fasta files for non-core gene families with >= 3 members
# For information on how to get monocopy core-gene families equivalent, 
# check the following directory in the same repository
# ../pangenomes-core_phylogenetic_trees/03_recombination_test.sh

# Prepare working directories
wDir="/my/path/to/working/directory"
panoctDir="$wDir/panoct_results"
nonCoreN3Dir="$panoctDir/nonCoreN3"
mkdir "$nonCoreN3Dir"

# Get entries n>=3:
cd "$nonCoreN3Dir"
cat $panoctDir/geneFamilies_panoct.txt | cut -f2 | sort | uniq > allLTs.list
cat $panoctDir/uniqueGenes_panoct.txt  | cut -f2 | sort | uniq > uniqlts.list
cat $panoctDir/coreGeneFamilies_panoct.txt  | cut -f2 | sort | uniq > corelts.list

# Get counts
wc -l *

# Get locus_tags that are neither from core genes nor unique genes and check counts
comm -23 allLTs.list uniqlts.list | sort | uniq > nonUniq.list.tmp 
comm -23 nonUniq.list.tmp corelts.list | sort | uniq > nonCore_nonUniqlts.list
wc -l nonCore_nonUniqlts.list 

# Get info for nonCore genes to get gene family and check counts
grep -Ff nonCore_nonUniqlts.list $panoctDir/geneFamilies_panoct.txt > nonCore_nonUniq.info
wc -l nonCore_nonUniq.info 

# Get gene families with >=3 entries
cat nonCore_nonUniq.info | cut -f1 | sort | uniq -c | perl -pe 's/\A +//g; s/ /\t/g;' | sort -nk1 | awk '$1 >= 3' | cut -f2 | sort | uniq > N3gf.list

# Map gene family ID to annotation; match by locus_tag. Then, check counts
# Matching a particular field in the second file (using '\t' delimiter and field 1 ):
awk -F '\t' 'FNR==NR {hash[$0]; next} $1 in hash' N3gf.list $panoctDir/geneFamilies_panoct.txt > nonCoreN3.info
cat nonCoreN3.info | cut -f1 | sort | uniq | wc -l 
wc -l N3gf.list 

# Delete intermediate files
rm *.list *tmp nonCore_nonUniq.info

# Extract fasta (fna) files
python split_faa_fna_ByPanocGeneFam.py -f $wDir/ -a $wDir/all_proteomes.faa -p $wDir/all_proteomes.attributes -g $nonCoreN3Dir/nonCoreN3.info -o $nonCoreN3Dir/


