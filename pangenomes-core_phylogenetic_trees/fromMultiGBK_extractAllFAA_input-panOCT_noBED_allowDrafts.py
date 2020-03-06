# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 16:01:24 2018

@author: mcastro
"""
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inGBK", help="COMPLETE path with input genbank file")
parser.add_argument("-g", "--genomeID", help="Genome identifier to be used; for example, strain name")
parser.add_argument("-o", "--outDir", help="COMPLETE path to output dir")
args = parser.parse_args()

#Extract CDS from GBK
def extractCDS(inGBK, genomeID, outDir):
    f=open(outDir+genomeID+".FAA","w")
    g=open(outDir+genomeID+".fasta","w")
    k=open(outDir+genomeID+".panoctInput","w")
    for record in SeqIO.parse(inGBK,"genbank"):
        accNumber=record.id
        g.write('>'+accNumber+"\n"+str(record.seq)+"\n")
        
        cdss = [feat.qualifiers for feat in record.features if feat.type == 'CDS']
        pos= [feat for feat in record.features if feat.type =='CDS']
        
        locus_pos={}
        
        exclude=[]
        
        for gen in pos:
            if 'pseudo' in gen.qualifiers.keys(): continue #skip pseudogenes
            elif 'ribosomal_slippage' in gen.qualifiers.keys(): continue #skip programmed frameshifts
            LT=gen.qualifiers['locus_tag'][0]
            if '>' in str(gen) or '<' in str(gen): continue #Omit incomplete genes
            elif 'location: join' in str(gen): #Join cases are special: Always start > end; and cannot be interchanged
                coords=str(gen)
                coords=coords.split('location')[1]
                coords=coords.split('qualifiers')[0]
                coords=coords.replace(': join{', '')
                coords=coords.replace('}\n', '')
                                
                a=coords.split('[')[1]
                a1=int(a.split(':')[0])
                a2=int(a.split(':')[1].split(']')[0])
                
                b=coords.split('[')[2]
                b1=int(b.split(':')[0])
                b2=int(b.split(':')[1].split(']')[0])
                
                #        [0, 775, 3823988, 3824908]
                end=sorted([a1,a2,b1,b2])[1] #775(a2) would be the end in this case 
                start=sorted([a1,a2,b1,b2])[2] #3823988(b1) would be the start                
                #end=a2
                #start=b1
                exclude.append(LT)
            else: #Normal cases
                start=int(gen.location.start)
                end=int(gen.location.end)
                
                if start > end:
                    tmp=start
                    start=end
                    end=tmp
            
            locus_pos.update({LT:[start, end]})
        
        for CDS in cdss:
            if 'locus_tag' in CDS and 'product' in CDS and 'translation' in CDS:
                LT=CDS['locus_tag'][0]
                if LT not in locus_pos: continue #skip pseudogenes and programmed frameshifts
                
                alias=LT
                f.write(">"+alias+" "+CDS['product'][0]+"\n"+CDS['translation'][0]+"\n")
                start=locus_pos[LT][0]
                end=locus_pos[LT][1]
                
                #Gene attribute file for panOCT
                # contig id, protein identifier (e.g. locus), 5'-coordinate, 3'-coordinate, annotation and genome identifier.
                k.write(accNumber+"\t"+alias+"\t"+str(start)+"\t"+str(end)+"\t"+CDS['product'][0]+"\t"+genomeID+"\n")
    
    f.close()
    g.close()
    k.close()
    return


extractCDS(args.inGBK, args.genomeID, args.outDir)



