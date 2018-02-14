#! /usr/bin/env python
import sys

'''
Get_Fasta_Seq.py [FASTA FILE] [SEQUENCE NAME]
Get_Fasta_Seq takes a FASTA file as its input, as well as the name of a specific sequence (e.g. chromosome)
Get_Fasta_Seq will print the sequence with the given name within the given FASTA file

Example: 
python Get_Fasta_Seq.py Yeast_Genome.fasta chrI
'''

file_path = sys.argv[1]
seq_name = sys.argv[2]

cur_seq_list = ['Sequence with the name {} not found!'.format(seq_name)]
cur_name = ''
f = open(file_path,'r')
for line in f:
    line = line.strip()
     
    if line.startswith('>'):
        cur_name = line[1:]
        if cur_name == seq_name:
            cur_seq_list = []    
    
    elif cur_name == seq_name:
        cur_seq_list.append(line)

print ''.join(cur_seq_list)
