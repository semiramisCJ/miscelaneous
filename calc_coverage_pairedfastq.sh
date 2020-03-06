# Calculate coverage for Illumina-only assemblies
# Case with paired-end reads in fastq format

# Our genome file is "AN3.fasta" and we are in the working directory
# Fastq read files are gzipped, but this is optional. 
# If files are not compressed, the commands stay the same; just provide the correct file extension.

#Index fasta file with bowtie2 
#Then, align reads and convert sam to bam with samtools
bowtie2-build AN3.fasta AN3.fasta
bowtie2 -x AN3.fasta -p 3 -1 AN3_R1.fastq.gz -2 AN3_R2.fastq.gz > AN3.sam #Align fastq gzipped files
samtools view -b -S AN3.sam > AN3.bam

#Then, concatenate bam files, sort them, and get coverage with bedtools
samtools sort AN3.bam > AN3.bam.sorted
bedtools genomecov -ibam AN3.bam.sorted -d -g AN3.fasta > AN3.depthCoverage.txt

#Print coverage to outfile
more AN3.depthCoverage.txt | awk '{sum+=$3;cnt++}END{print str" "sum/cnt" "sum}' str="AN3" > AN3.coverageTable.txt

