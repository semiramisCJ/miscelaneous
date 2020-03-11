##################### codeml for coreProts ###################
gf=$1
wDirBase="/my/path/to/working/directory"
alnDirBase="$wDirBase/panoct_results"
templateDir="/my/path/to/working/directory/templatesCodeML/"
currCase="coreProts"
wDir="$alnDirBase/$currCase"
phyDir="$wDir/phyDir"
outDir="$wDir/codeml_selection"

#Run analysis for each gene family
cd "$wDir" #Go to working dir
models=( M1 M2 ) 

#Convert aln to paml-compatible phylip format
perl /my/installation/path/Fasta2Phylip.pl $gf.fas $gf.phy.tmp
perl -pe 's/\t/\n/g' $gf.phy.tmp > $gf.phy

#Convert tree to newick-compatible format 
perl -pe 's/\)\d+/\)/g;' $phyDir/$gf.fas.treefile > $gf.tree

seqFileName="$gf.phy"
treeFileName="$gf.tree"

#Create empty table to write lnL
> $gf.lnL.tab

#Run models	
for model in "${models[@]}"; do
	resultsFileName="$gf.$model.results.txt"
	#Edit model control file
	cat $templateDir/$model.codeml.ctl | perl -spe 's/seqfile\.txt/$newSeq/; s/treefile\.txt/$newTree/; s/results\.txt/$newRes/' -- -newSeq=$seqFileName -newTree=$treeFileName -newRes=$resultsFileName > codeml.ctl

	#Run codeml
	codeml

	#Delete old control file and temporal files
	rm codeml.ctl
	rm $gf.phy.tmp

	#Get lnL from each model for each gene family and write it to table
	lnL=$(grep lnL $gf.$model.results.txt | perl -pe 's/\):/\t/g; s/\+/\t+/g;' | cut -f2 | perl -pe 's/ +//g;')
	echo -e "$model\t$lnL" >> $gf.lnL.tab
done

##################################################

