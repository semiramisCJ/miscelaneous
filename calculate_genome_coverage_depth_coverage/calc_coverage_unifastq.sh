# Calculate coverage for hybrid assemblies
# Case with paired-end reads in fastq format and long reads also in fastq format (as with PacBio RSII filtered_reads downloaded from SMRT_portal)

# Our genome file is "AN54.fasta" and we are in the working directory
# Fastq read files are gzipped, but this is optional. 
# If files are compressed, the commands stay the same; just provide the correct file extension.

# First of all, concatenate all fastq reads into one big fastq file
cat *.fastq > ALL_reads.fastq

# Then, index fasta file with bwa, align reads and convert sam to bam with samtools
bwa index AN54.fasta
bwa mem -t 6 AN54.fasta ALL_reads.fastq > gAN54.aln.sam
samtools view -b -S gAN54.aln.sam > gAN54.aln.bam

# Now, sort bam file and get coverage with bedtools
samtools sort gAN54.aln.bam > gAN54.aln.bam.sorted
bedtools genomecov -ibam gAN54.aln.bam.sorted -d -g AN54.fasta > gAN54.depthCoverage.txt

# Finally, print coverage to outfile
more gAN54.depthCoverage.txt | awk '{sum+=$3;cnt++}END{print str" "sum/cnt" "sum}' str="AN54" > gAN54.coverageTable.txt


