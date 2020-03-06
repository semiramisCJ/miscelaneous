# Get homologous groups with get_homologues pipeline and use get_phylomarkers

# Depending on the installation, we might need to activate a specific conda environment to avoid compatibility issues
source activate get_phylomarkers

# get_homologues pipeline needs genbank files ending with *.gbk, so keep this in mind when preparing the input
# genbank files should have both DNA sequence and annotation

# Create working directories, copy and rename files
wDirBase="/my/path/to/working/directory/"
genomeDir="genomes_Ahaem"
cd "$wDirBase"
cp /my/path/to/annotation/*.gbff .
rename gbff gbk *.gbff
mkdir "$genomeDir" # This is because we need a separate directory for the input genbanks
mv *.gbk $genomeDir/


# Run get_homologues with different algorithms for further comparison. We decided to keep the intersection of all methods
#Global params
#-C min %coverage in BLAST pairwise alignments
#-S min %sequence identity in BLAST query/subj pairs

# Bi-directional best hit
get_homologues.pl -d genomes_Ahaem -C 80 -S 90 &> log.BDBH_def 

# PFAM
get_homologues.pl -d genomes_Ahaem -D -C 80 -S 90 &> log.pfam

# COG triangles
get_homologues.pl -d genomes_Ahaem -G -t 0 -c -C 80 -S 90 &> log.cog 

# OrthoMCL
get_homologues.pl -d genomes_Ahaem -M -t 0 -c -C 80 -S 90 &> log.orthomcl


# After all runs are complete, list directories to manually give them to the intersection comparison program (compare_clusters.pl)
wDirBase="/my/path/to/working/directory/"
homologuesDir="$wDirBase/genomes_Ahaem_homologues" #Its name is $genomeDir + "_homologues"
cd "$homologuesDir"
find . -type d
.
./tmp
./sz1652_f0_alltaxa_algBDBH_e0_C80_S90_
./sz1652_f0_alltaxa_algBDBH_Pfam_e0_C80_S90_
./sz1652_f0_0taxa_algCOG_e0_C80_S90_
./sz1652_f0_0taxa_algOMCL_e0_C80_S90_

# Compare clusters for pangenome; first use proteins; then use nucleotides. The pangenome matrix is one of its main outputs
compare_clusters.pl -d sz1652_f0_0taxa_algCOG_e0_C80_S90_,sz1652_f0_0taxa_algOMCL_e0_C80_S90_ -o intersect_pan -t 0
compare_clusters.pl -d sz1652_f0_0taxa_algCOG_e0_C80_S90_,sz1652_f0_0taxa_algOMCL_e0_C80_S90_ -o intersect_pan -t 0 -n -m

# Compare clusters for core genome; get_phylomarkers needs this output
compare_clusters.pl -d sz1652_f0_alltaxa_algBDBH_e0_C80_S90_,sz1652_f0_alltaxa_algBDBH_Pfam_e0_C80_S90_,sz1652_f0_0taxa_algCOG_e0_C80_S90_,sz1652_f0_0taxa_algOMCL_e0_C80_S90_ -o intersect_core -t 12 
compare_clusters.pl -d sz1652_f0_alltaxa_algBDBH_e0_C80_S90_,sz1652_f0_alltaxa_algBDBH_Pfam_e0_C80_S90_,sz1652_f0_0taxa_algCOG_e0_C80_S90_,sz1652_f0_0taxa_algOMCL_e0_C80_S90_ -o intersect_core -t 12 -n -m


# Run get_phylomarkers on Phylogenomics mode; for example to use the gene trees for species phylogenies
cd $homologuesDir/intersect_core
run_get_phylomarkers_pipeline.sh -R 1 -t DNA

# Run get_phylomarkers on Population genetics mode; for intra-specific comparisons. It outputs a table with Tajima's D and other useful metrics
cd $homologuesDir/intersect_core
run_get_phylomarkers_pipeline.sh -R 2 -t DNA

# If we only want to get pangenome plots, we need to re-run get_homologues pipeline with flag -c (genome compositional analysis)
# The flag -c will create the input files for plot_pancore_matrix.pl, otherwise plot_pancore_matrix will not run.
cd "$wDirBase"
get_homologues.pl -d $genomeDir -c
cd $homologuesDir
plot_pancore_matrix.pl -i pan_genome_algBDBH.tab -f pan
plot_pancore_matrix.pl -i core_genome_algBDBH.tab -f core_both


# Leave the environment
conda deactivate


