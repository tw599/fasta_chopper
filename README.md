# fasta_chopper
Python script for 'chopping up' DNA sequences/reads/regions into uniform-sized fragments. This is practically useful for oligo tiling designs that are often used in large scale oligo synthesis and assembly experiements. All sequences generated contain fixed overlap regions for homology cloning/contig assembly etc. 

# Hardware/Software requirements

- x86-64 compatible processors
- 64 bit Linux or Mac OS X

# Usage 

Current version needs two arguments for execution.
"input.fasta" and "output.fasta" filenames & path need to be entered for successful execution of the script.

Script for execution of script from terminal would be as follows:

```
python fasta_chopper.py input.fasta output.fasta
```

Now includes options to enter the following parameters:

*barcode_file* can be included to add in a unique barcode sequence at the 5' end of each tiled oligo

*site1, site2* : allow inclusion of unique primer sites for amplification and PCR purposes

*--window_size* parameter can now be specified for determining the length of the tiled oligos

*--step_size* parameter can now be specified to determine the overlap regions between sequentially 'chopped' oligos

# Mandatory Parameters

*input.fasta* and *output.fasta* are the only two mandatory requirements

# Author
Tushar Warrier 
