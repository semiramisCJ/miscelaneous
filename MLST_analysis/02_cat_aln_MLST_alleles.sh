# Align each locus, filter gaps, and concatenate alignments

# First, create fasta files for each profile in MLST DB
working_dir="/my/path/to/working/directory"
db_seq_dir="$working_dir/pubmlst_cron_sch1/"
output_dir="$working_dir/pubmlst_cron_sch1/alleles/"
profile_file_name="profile_scheme1.tab"
mkdir "$output_dir"
python build_MLST_profile_fasta_for_aln.py --db_dir $db_seq_dir --profile_file_name $profile_file_name --output_dir $output_dir

# Check results
cd "$output_dir"
ls 

# Now, format query results; delete extra strings and new lines in fasta
query_seq_dir="$working_dir/alleles_cronobacter_1"
cd "$query_seq_dir"
ls *.fasta > locus.list
while read locus
do
	perl -pe ' s/\n//; s/>/\n>/g; s/>.+;(.+);.+/>$1\n/g;' $locus > parsed.$locus
done < locus.list

# Then, concatenate query and DB sequences
db_seq_dir="$working_dir/pubmlst_cron_sch1/alleles"
query_seq_dir="$working_dir/alleles_cronobacter_1"
aln_dir="$working_dir/MLST_sequences_for_aln"
mkdir "$aln_dir"
cd "$aln_dir"
while read locus
do
	cat $db_seq_dir/$locus $query_seq_dir/parsed.$locus > $aln_dir/$locus
done < $query_seq_dir/locus.list


# Rename files and create loci list 
cd "$aln_dir"
rename fasta fna *.fasta
ls *.fna | perl -pe 's/\.fna//g;' > locus.list

# Create directories for aligned sequences: a final directory and a temporal directory
mkdir alignments_ungapped alignments_tmp

### G-INS-i assumes that entire region can be aligned and tries to align them globally using the Needleman-Wunsch algorithm; that is, a set of sequences of one domain must be extracted by truncating flanking sequences
while read locus
do
	# Check orientation of sequences: all in forward
	mafft --adjustdirection $locus.fna > alignments_tmp/tmp 
	
	# Delete extra letters on flipped sequences
	perl -pe 's/>_R_/>/g;' alignments_tmp/tmp > alignments_tmp/oriented.$locus 
	
	# Align with G-INS-i algorithm
	mafft --globalpair --maxiterate 1000 --thread 5 --op 2 alignments_tmp/oriented.$locus > alignments_tmp/$locus.fas 
	#--op # :         Gap opening penalty, default: 1.53
	
	# Delete gapped columns with trimal
	/my/installation/path/trimal -in alignments_tmp/$locus.fas -out alignments_ungapped/$locus.ungapped.fas -nogaps 
	#-nogaps                  Remove all positions with gaps in the alignment.
done < locus.list


# Cat alignments with FASTconcat
cd $aln_dir/alignments_ungapped
perl /my/installation/path/FASconCAT-G.pl -s -n -p -p ##And press enter
#-f Defined Input files [All Input Files]
#-s Concatenation
#-n Print NEXUS Output (Block)
# -p -p Print PHYLIP Output (Relaxed)

#Rename files
rename FcC_supermatrix all_mlst_locus FcC_supermatrix.*



