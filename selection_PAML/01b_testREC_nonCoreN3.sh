# Align CDS guided by proteins & unify names & perform REC test

# Go to working directory
wDir="/my/path/to/working/directory"
panoctDir="$wDir/panoct_results"
nonCoreN3Dir="$panoctDir/nonCoreN3"
recDir="$nonCoreN3Dir/withRECsignal/"
mkdir $recDir
cd "$nonCoreN3Dir"

# Create empty files for REC tests
> unableToTest.txt
> significative.phiResults.all
cut -f1 nonCoreN3.info | sort | uniq > groups.list
while read group
do
	# Check if file exists
	if [ ! -f "$group".faa ]; then
	    echo "File not found for group $group"
	    continue #Skip
	fi
	
	# Align proteins with clustal omega
	clustalo -i $group.faa -t Protein -o aln_$group.fasta
		
	# Run revtrans
	python /my/installation/path/revtrans.py $group.fna aln_$group.fasta -match name -O fasta RevTransAln_$group.fa
	# Delete protein aln, proteins and genes 
	rm aln_$group.fasta $group.faa $group.fna
	
	# Unify names & delete orig file
	perl -pe 's/(>.+) .+/$1/; s/\.gbff//g; s/\.gbf//g;' RevTransAln_$group.fa > $group.fas
	rm RevTransAln_$group.fa

	# Perform REC test with phi-pack
	/my/installation/path/Phi -f $group.fas -p -w 50 > results.tmp	
	#  -f: Filename = FASTA format
	#  -p: [#] = PHI permutation test [default = FALSE, #=1000]
	#  -w: # = Change default window size [default w = 100]
	
	# Test if the exit code of phi phack ($?) is 1; that means the program couldn't run; Write to outfile a list of those problematic alns
	if [ "$?" == "1" ] ; then echo $group >> unableToTest.txt; fi
	
	# Get and reformat custom columns
	cat results.tmp | grep ':' | grep -v Writing | grep -v Done | grep -v Estimated | perl -pe 's/: +/\t/g; s/PHI +\(/PHI\(/g; s/ +\(/\t/g; s/\)\n/\n/g; s/ +//g;' | awk -v name="$group" -F $'\t' '{ print name,"\t",$1,"\t",$2,"\t",$3 }' | cut -f1-3 > $group.phiResults
	
	# Print significative results
	awk ' $3 < 0.05 ' $group.phiResults | grep -v "\-\-" >> significative.phiResults.all
	
	# Delete intermediate files	
	rm $group.phiResults results.tmp

done < groups.list

# Explore final output
cat significative.phiResults.all
grep 'Permutation' significative.phiResults.all | cut -f1 | sort | uniq > withRECsignal.list

# Move alns with REC signal to other directory
while read group
do
	mv $group.fas $recDir/
done < withRECsignal.list

