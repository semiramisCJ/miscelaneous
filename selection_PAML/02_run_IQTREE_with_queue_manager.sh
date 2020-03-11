# Use IQ-TREE suite to get best model and run phylogeny for each gene family
# Use the queue manager available in the server

# Go to working directory
wDir="/my/path/to/working/directory"
panoctDir="$wDir/panoct_results"
nonCoreN3Dir="$panoctDir/nonCoreN3"
coreDir="$panoctDir/coreProts"

# List aln_fasta files without rec signals for core and non-core gene families
# Then, print task files

# For core gene families
cd "$coreDir"
ls *.fas > noREC.coreProts.list
ls *.fas | perl -pe 's/\.fas//g;' > noREC.coreProts.gf.list
while read currGF
do
	echo "bash coreProts_IQtree.sh $currGF" >> taskList.txt
done < noREC.coreProts.gf.list

# For non-core gene families
cd "$nonCoreN3Dir"
ls *.fas > noREC.nonCoreN3.list
ls *.fas | perl -pe 's/\.fas//g;' > noREC.nonCoreN3.gf.list
while read currGF
do
	echo "bash nonCoreN3_IQtree.sh $currGF" >> taskList.txt
done < noREC.nonCoreN3.gf.list

# Before submitting array jobs,
# edit jdl files and bash scripts
# The scripts are in the folder "files_for_queue_manager"
# To obtain the queue name you will work with, contact your system administrator
# Then, check the number of entries to define batch size
# The following case is just an example

# For core gene families
cd "$coreDir"
wc -l noREC.coreProts.gf.list 
qsub -t 1-10 iqTree_coreProts.jdl
qsub -t 11-100 iqTree_coreProts.jdl
qsub -t 101-200 iqTree_coreProts.jdl
qsub -t 201-300 iqTree_coreProts.jdl
qsub -t 301-400 iqTree_coreProts.jdl
qsub -t 401-500 iqTree_coreProts.jdl
qsub -t 501-600 iqTree_coreProts.jdl
qsub -t 601-700 iqTree_coreProts.jdl
qsub -t 701-800 iqTree_coreProts.jdl
qsub -t 801-900 iqTree_coreProts.jdl
qsub -t 901-921 iqTree_coreProts.jdl

# For non-core gene families
cd "$nonCoreN3Dir"
wc -l noREC.nonCoreN3.gf.list 
qsub -t 1-500 iqTree_nonCoreN3.jdl
qsub -t 501-1000 iqTree_nonCoreN3.jdl
qsub -t 1001-1500 iqTree_nonCoreN3.jdl
qsub -t 1501-1750 iqTree_nonCoreN3.jdl
qsub -t 1751-2000 iqTree_nonCoreN3.jdl
qsub -t 2001-2500 iqTree_nonCoreN3.jdl
qsub -t 2501-3000 iqTree_nonCoreN3.jdl
qsub -t 3001-3401 iqTree_nonCoreN3.jdl


# * Note: monitor jobs with 
qstat
qstat -u \*


