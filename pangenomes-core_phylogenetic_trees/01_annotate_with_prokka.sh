#Annotate all genomes with Prokka

wDirBase="/my/path/to/working/directory/"
cd "$wDirBase"

ls *.fasta > genomes.list
while read fname
do
	prokka --cpus 8 --kingdom Bacteria --gram neg --genus Acinetobacter --strain $fname --usegenus --addgenes --locustag $fname --prefix $fname --outdir $fname.prokka $fname
done < genomes.list

