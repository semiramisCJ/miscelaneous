# miscelaneous scripts
Diverse scripts for comparative genomics on bacteria.

Example file names and paths are provided; if you use these scripts, please make the appropriate changes. Also, check if the required software is available in your computer, and if the number of threads used make sense for your particular situation.

The examples are done for different genomes / genus / species

Name consistency is inside folders; variable names are used for most of the parts for easy editing.


# Dependencies
Third-party software is called under the "/my/installation/path/" or directly; depends on the installation we had at the moment.

A (non-exhaustive) list of **dependencies that vary by script**:
* bowtie2
* bwa
* samtools
* bamtools
* blast+
* python3 with argrparse, BioPython, networkx, click
* perl
* prokka
* panoct
* RevTrans
* Clustal Omega
* PhiPack
* FASconCAT-G
* RAxML
* GET_HOMOLOGUES
* GET_PHYLOMARKERS
* MAFFT
* IQ-TREE
* Trimal 
* R with MLSTaR, purrr, readr, dplyr
* Fasta2Phylip.pl obtained from https://github.com/josephhughes/Sequence-manipulation/blob/master/Fasta2Phylip.pl


# License
MIT License

Copyright (c) [2020] [Semiramis Castro-Jaimes]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
