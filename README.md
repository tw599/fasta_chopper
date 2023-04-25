# fasta_chopper
Python script for 'chopping up' DNA sequences/reads/regions into uniform-sized fragments. This is practically useful for oligo tiling designs that are often used in large scale oligo synthesis and assembly experiements. All sequences generated contain fixed overlap regions for homology cloning/contig assembly etc. 

Author:
Tushar Warrier 
tw599@cam.ac.uk
https://github.com/tw599/fasta_chopper/issues

Hardware/Software requirements:
x86-64 compatible processors
64 bit Linux or Mac OS X

MANUAL:

The script is very easy to use. Current version needs two arguments for executtion.
"input.fasta" and "output.fasta" filenames need to be entered for successful execution of the script.

Script for execution of script from terminal would be as follows:

python fasta_chopper.py input.fasta output.fasta

The current default length of chopped fragments is 220 bp with an overlap of 20 bp. The last chopped fragment for each record of the FASTA file (each DNA sequence) is returned as a 220-bp sequence which has a >20 bp overlap sequence.
This only occurs in the event that:

length of last chopped sequence/overhang < 220 bp (chopping length)


For editing chop length:

Edit lines 11 and 12 of fasta_chopper.py file:

window_size = 220 # length of each sequence
step_size = 200 # distance between the start positions of adjacent sequences

NOTE: Subsequent versions will allow for entry of the window_size and step_size as user-determined parameters in the command-line/terminal.
NOTE: Would recommend users to add a unique classifier to each DNA sequence (record) in their fasta file. This will allow for easier reading of the output DNA chunks.

