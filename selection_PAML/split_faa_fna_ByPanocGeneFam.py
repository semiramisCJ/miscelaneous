#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 13:43:04 2019

@author: mcastro
"""
from Bio import SeqIO
from Bio import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--faDir", help="COMPLETE path to the input directory where *.fasta files are located")
parser.add_argument("-a", "--allProteomesFile", help="COMPLETE path to the faa multifasta file used as input for blast and panoct")
parser.add_argument("-p", "--panoctAttributesFile", help="COMPLETE path to the attributes file used as input for panoct")
parser.add_argument("-g", "--geneFamFile", help="COMPLETE path to parsed file with gene family + LT ")
parser.add_argument("-o", "--outDir", help="COMPLETE path to the outDir")
args = parser.parse_args()



def getColumns( infileName ):
    "Parses a tab-delimited file and stores it in a list of lists. This allows to easily access the columns"
    data= open(infileName,"rU")
    lines= data.readlines()
    precolumns= [l.strip("\n") for l in lines]
    data.close()
    columns= [l.split("\t") for l in precolumns]
    return columns



def splitFams_FAA_FNA(faDir, allProteomesFile, panoctAttributesFile, geneFamFile, outDir):
    #Load proteomes
    protDict=SeqIO.to_dict(SeqIO.parse(allProteomesFile, "fasta"))
    
    #Read start-end coords, contig and strain from panoct input
    allAttributes=getColumns(panoctAttributesFile)
    #LT:[strain, start, end, contig]
    allAttributes={c[1]:[c[5], c[2], c[3], c[0]] for c in allAttributes}
    
    
    #Load fna file to extract nts and print CDS of the same gene family to the same outfile
    famDict={}
    families=getColumns(geneFamFile)
    for line in families:
        fam=line[0]
        LT=line[1]
                
        if fam not in famDict:
            famDict.update({fam:{}})
        
        famDict[fam].update({LT:allAttributes[LT]})
    
    
    #Iterate through each family with len > 1
    for fam in famDict:
        if len(famDict[fam].keys()) == 1: continue
        fna=open(outDir+str(fam)+".fna", "w")
        faa=open(outDir+str(fam)+".faa", "w")
        for LT in famDict[fam]:
            strain=famDict[fam][LT][0]
            start=int(famDict[fam][LT][1])
            end=int(famDict[fam][LT][2])
            contig=famDict[fam][LT][3]
            
            #Find the locus_tag in the protDir and write it to outfile
            faa.write('>'+strain+" "+LT+"\n"+str(protDict[LT].seq)+"\n")
            
            #Open the strain's fasta file to find the desired contig
            genome=SeqIO.parse(faDir+strain+".fasta", "fasta")
            for record in genome:
                if record.id == contig:
                    if start < end:
                        gene=str(record.seq)[start:end]
                    else:
                        part1=str(record.seq)[start:]
                        part2=str(record.seq)[:end]
                        gene=part1+part2
                    break
            
            #Get sequences in the same orientation (all in forward)
            gene=Seq.Seq(gene)
            try:
                gene.translate(table="Bacterial", cds=True)
            except:
                gene=gene.reverse_complement()
            
            gene=str(gene)
            fna.write('>'+strain+" "+LT+"\n"+gene+"\n")
        
        faa.close()
        fna.close()
    
    return
    
splitFams_FAA_FNA(args.faDir, args.allProteomesFile, args.panoctAttributesFile, 
                  args.geneFamFile, args.outDir) 
 
