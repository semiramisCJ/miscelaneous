#15 Mar 2018
#Calculate coverage for hybrid assemblies. Case with all reads in fastq format (Illumina + PacBio RSII)

###Run this in chichen
cd /space21/PGE/mcastro/ensamble_AN54

#Our genome file is "AN54_bestDraft_Mar2018.fasta"

#Cat all fastq reads into one big fastq file
cat *.fastq > ALL_reads.fastq

#Index fasta file with bwa, convert sam to bam with samtools, sort it, and get coverage with bedtools
bwa index AN54_bestDraft_Mar2018.fasta
bwa mem -t 6 AN54_bestDraft_Mar2018.fasta ALL_reads.fastq > gAN54.aln.sam
samtools view -b -S gAN54.aln.sam > gAN54.aln.bam
samtools sort gAN54.aln.bam > gAN54.aln.bam.sorted
bedtools genomecov -ibam gAN54.aln.bam.sorted -d -g AN54_bestDraft_Mar2018.fasta > gAN54.depthCoverage.txt

#Print coverage to outfile
more gAN54.depthCoverage.txt | awk '{sum+=$3;cnt++}END{print str" "sum/cnt" "sum}' str="gAN54" > gAN54.coverageTable.txt


