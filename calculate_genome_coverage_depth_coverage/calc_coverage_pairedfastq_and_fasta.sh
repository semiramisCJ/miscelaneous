# Calculate coverage for hybrid assemblies
# Case with paired-end reads in fastq format and long reads in fasta format (as with PacBio Sequel)

# Our genome file is "AN3.fasta" and we are in the working directory
# Fastq read files are gzipped, but this is optional. 
# If files are not compressed, the commands stay the same; just provide the correct file extension.

# First, index fasta file with bowtie2 
# Then, align reads and convert sam to bam with samtools
bowtie2-build AN3.fasta AN3.fasta
bowtie2 -x AN3.fasta -p 3 -1 AN3_R1.fastq.gz -2 AN3_R2.fastq.gz > AN3_illum.sam # Align fastq gzipped files
bowtie2 -x AN3.fasta -p 3 -f -U AN3_Sequel.fasta.gz > AN3_sequel.sam # Align fasta gzipped files
samtools view -b -S AN3_illum.sam > AN3_illum.bam
samtools view -b -S AN3_sequel.sam > AN3_sequel.bam

# Now, concatenate bam files, sort them, and get coverage with bedtools
samtools cat AN3_illum.bam AN3_sequel.bam > AN3.bam
samtools sort AN3.bam > AN3.bam.sorted
bedtools genomecov -ibam AN3.bam.sorted -d -g AN3.fasta > AN3.depthCoverage.txt

# Finally, print coverage to outfile
more AN3.depthCoverage.txt | awk '{sum+=$3;cnt++}END{print str" "sum/cnt" "sum}' str="AN3" > AN3.coverageTable.txt

