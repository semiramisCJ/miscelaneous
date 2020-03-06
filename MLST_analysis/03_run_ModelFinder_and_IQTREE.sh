# Search best model and run tree search with ModelFinder and IQ-TREE

# Go to directory with ungapped, concatenated alignment
working_dir="/my/path/to/working/directory"
aln_dir="$working_dir/MLST_sequences_for_aln"
cd $aln_dir/alignments_ungapped

# Run ModelFinder
nohup /my/installation/path/iqtree -s all_mlst_locus.fas -m MFP &

# See ModelFinder results; the result will be passed to the next command with -m
grep 'Model of substitution' all_mlst_locus.fas.iqtree

# Run tree reconstruction with selected model with iqtree
nohup /my/installation/path/iqtree -s all_mlst_locus.fas -m GTR+F+I+G4 -b 100 -nt 3 -redo &
#nt number of threads

# Check progress with grep
grep 'BOOTSTRAP REPLICATE NUMBER' all_mlst_locus.fas.log


