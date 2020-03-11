# Run selection tests with PAML suite
# *Note Fasta2Phylip.pl was obtained from https://github.com/josephhughes/Sequence-manipulation/blob/master/Fasta2Phylip.pl

# Go to working directory
wDir="/my/path/to/working/directory"
panoctDir="$wDir/panoct_results"
nonCoreN3Dir="$panoctDir/nonCoreN3"
coreDir="$panoctDir/coreProts"

# Remove previous task files; then, print new task files

# For core gene families
cd "$coreDir"
rm taskList.txt
> taskList.txt
while read currGF
do
	echo "bash coreProts_runCodeml.sh $currGF" >> taskList.txt
done < noREC.coreProts.gf.list

# For non-core gene families
cd "$nonCoreN3Dir"
rm taskList.txt
> taskList.txt
#ls *.fas | perl -pe 's/\.fas//g;' > noREC.nonCoreN3.gf.list
while read currGF
do
	echo "bash nonCoreN3_runCodeml.sh $currGF" >> taskList.txt
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
qsub -t 1-100 codeml_coreProts.jdl
qsub -t 83-100 codeml_coreProts.jdl
qsub -t 101-200 codeml_coreProts.jdl
qsub -t 201-300 codeml_coreProts.jdl
qsub -t 301-400 codeml_coreProts.jdl
qsub -t 401-500 codeml_coreProts.jdl
qsub -t 501-600 codeml_coreProts.jdl
qsub -t 601-700 codeml_coreProts.jdl
qsub -t 701-800 codeml_coreProts.jdl
qsub -t 801-900 codeml_coreProts.jdl
qsub -t 901-921 codeml_coreProts.jdl


# For non-core gene families
cd "$nonCoreN3Dir"
wc -l noREC.nonCoreN3.gf.list 
qsub -t 1-500 codeml_nonCoreN3.jdl
qsub -t 501-1000 codeml_nonCoreN3.jdl
qsub -t 555-1000 codeml_nonCoreN3.jdl
qsub -t 1001-1500 codeml_nonCoreN3.jdl
qsub -t 1501-1750 codeml_nonCoreN3.jdl
qsub -t 1751-2000 codeml_nonCoreN3.jdl
qsub -t 2001-2500 codeml_nonCoreN3.jdl
qsub -t 2501-3000 codeml_nonCoreN3.jdl
qsub -t 3001-3401 codeml_nonCoreN3.jdl


# * Note: monitor jobs with 
qstat
qstat -u \*

