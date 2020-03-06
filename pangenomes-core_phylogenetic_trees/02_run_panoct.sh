# Create panoct input, run blast search, and run panoct


# Create panoct working directory
wDirBase="/my/path/to/working/directory"
wDir="$wDirBase/panoct_input/"
mkdir "$wDir"

# Copy gbff files to new working dir
cd "$wDir"
cp $wDirBase/*.prokka/*.gbk .
rename gbk gbf *gbk

# Parse gbff files to obtain faa, fasta, and panoct input files from complete proteomes
ls *.gbf > gbkList.txt
while read currGBK
do
	python fromMultiGBK_extractAllFAA_input-panOCT_noBED_allowDrafts.py -i $wDir/$currGBK -g $currGBK -o $wDir
done < gbkList.txt

# Cat all the FAA files and delete intermediate files
cat *.FAA > all_proteomes.faa
rm *.FAA 

# Run blastp all vs all
makeblastdb -in all_proteomes.faa -dbtype prot
blastp -db all_proteomes.faa -query all_proteomes.faa -num_threads 10 -evalue 0.0001 -outfmt 6 -out all_vs_all.blastp.tab

# Run panoct
cat *.panoctInput > all_proteomes.attributes #Cat the genome attribute files
cut -f6 all_proteomes.attributes | sort | uniq > genomesList_panoct.txt #Create genomes list file
panoctDir="$wDir/panoct_results"
mkdir "$panoctDir"
cd "$panoctDir"
perl /my/installation/path/panoct.pl -t ../all_vs_all.blastp.tab -f ../genomesList_panoct.txt -g ../all_proteomes.attributes -P ../all_proteomes.faa -S Y -L 0.80 -i 65 -V Y -T -c 100

# Parse panoct results to get coreOrth and single genes
python parsePanOCTresults-allowDrafts_addcoreGFlist.py -i $wDir/ -p $panoctDir/ -s $wDir/genomesList_panoct.txt

# Extract genes by gene family; only for core genes
coreDir="$panoctDir/coreProts"
mkdir "$coreDir"
python split_faa_fna_ByPanocGeneFam.py -f $wDir/ -a $wDir/all_proteomes.faa -p $wDir/all_proteomes.attributes -g $panoctDir/coreGeneFamilies_panoct.txt -o $coreDir/

