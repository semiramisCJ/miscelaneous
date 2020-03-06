#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:10:02 2020

@author: semiramis
"""

import click
from Bio import SeqIO

def get_columns( infileName ):
    """Parses a tab-delimited file and stores it in a list of lists. 
    This allows to easily access the columns"""
    
    data = open(infileName,"r")
    lines = data.readlines()
    precolumns = [l.strip("\n") for l in lines]
    data.close()
    columns = [l.split("\t") for l in precolumns]
    return columns

@click.command()
@click.option("--db_dir", help="Complete path input directory with separate fasta alleles and table of profiles")
@click.option("--profile_file_name", help="File name of profile table")
@click.option("--output_dir", help="Complete path to writable output directory")
def build_fasta_from_profiles(db_dir, profile_file_name, output_dir):
    """Main function"""
    
    # Load profile table and get allele names
    profiles = get_columns(db_dir+profile_file_name)
    loci_list = profiles[0][1:-1]
    
    # Create a dictionary of profiles
    profiles = profiles[1:]
    profiles = {c[0]: c[1:-1] for c in profiles if c[0] != '' }
    
    # Build a dictionary with all allele entries for all loci in scheme
    # And create an empty fasta file for each locus in output dir
    seq_dict={}
    for locus in loci_list:
        f=open(output_dir+locus+".fasta", "w")
        f.close()
        locus_fasta = SeqIO.parse(db_dir+locus+".fas", "fasta")
        for allele in locus_fasta:
            seq_dict.update({allele.id: allele.seq})
    
    # Write sequences for each for each ST for each locus
    for ST in profiles:
        for i in range(len(loci_list)):
            locus = loci_list[i]
            allele = profiles[ST][i]
            
            with open(output_dir+locus+".fasta", "a") as f:
                f.write(">ST" + ST + "\n")
                f.write(str(seq_dict[locus+"_"+allele]) + "\n")
    
    return

# Call function
build_fasta_from_profiles()
