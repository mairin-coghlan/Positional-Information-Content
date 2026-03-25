import pandas as pd
import numpy as np
import os
from Bio import motifs
from Bio.Seq import Seq
from math import log2

file = '14test.txt'
os.chdir('/stratagem/processed_16s/')

#processing text file into list
def text_to_list(file):
    sequences = []

    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                sequences.append(line)
    
    return sequences

#print(*text_to_list(file), sep = '\n')


#turning barcodes into Seq objects to be used in motif creation (needed for bio.motif functions)
def process_barcodes(barcodes):
    processed_barcodes = []

    for barcode in barcodes:
        i = Seq(barcode)
        processed_barcodes.append(i)
    return processed_barcodes

#calculating positional info content using given list of motifs
def calculate_pic(motif):

    position_ic = []

    for i in range(motif.length): #for number of positions across barcodes of same length

        distinct_bases = 0 #initializing count

        for b in ['A', 'T', 'C', 'G']:

            if motif.counts[b][i] > 0: #if the base appears at least once at this position, increase count of distinct bases
                distinct_bases += 1

            if distinct_bases > 1: #if there is more than one distinct base at this position, calculate the information content using log2 of the number of distinct bases
                p_ic_val = log2(distinct_bases)

            else:
                p_ic_val = 0

        position_ic.append('pos ' + str(i+1) + ': ' + str(p_ic_val) + ' bits')
    
    return position_ic


#calling script
motif = motifs.create(process_barcodes(text_to_list(file)), alphabet = 'ATCG')
print(*calculate_pic(motif), sep = '\n')

#QC case
#barcodes = ['AAA', 'ACA', 'AGA', 'ATA', 'AAC', 'ACC', 'AGC', 'ATC', 'AAG', 'ACG', 'AGG', 'ATG', 'AAC', 'ACC', 'AGC', 'ATC']

#motif = motifs.create(process_barcodes(barcodes), alphabet='ATCG')

#print(motif.counts) #outputs matrix of how many times each base appears at each position
#print(motif.length) 
#print(motif.counts['A']) #outputs list of how many times A appears at each position

#print(*calculate_pic(motif), sep='\n')



