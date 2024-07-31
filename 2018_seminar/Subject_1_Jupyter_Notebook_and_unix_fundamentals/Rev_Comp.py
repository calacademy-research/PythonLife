#!/usr/bin/env python
import sys
'''
Rev_Comp.py [DNA sequence]
Rev_Comp.py takes a DNA sequence as its input
It will print the reverse complement of the sequence to the terminal

If no argument is given, Rev_Comp.py will use stdin (standard in, aka the output from a pipe)

Examples:
python Rev_Comp.py ATAGAG
cat file.seq | python Rev_Comp.py
'''


if len(sys.argv) == 1:
    input_string = sys.stdin.read()
else:
    input_string = sys.argv[1].strip()

string_list = []

for char in input_string:
    if char == 'A':
        string_list.append('T')
    elif char == 'C':
        string_list.append('G')
    elif char == 'G':
        string_list.append('C')
    elif char == 'T':
        string_list.append('A')    
    elif char == 'a':
        string_list.append('t')
    elif char == 'c':
        string_list.append('g')
    elif char == 'g':
        string_list.append('c')
    elif char == 't':
        string_list.append('a') 
    elif char == '' or char == '\n':
        string_list.append('') 
        break
    else:
        string_list.append('N')

print ''.join(string_list[::-1])
