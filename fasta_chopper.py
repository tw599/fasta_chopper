import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Read input and output files
input_file = sys.argv[1]
output_file = sys.argv[2]

# Define parameters
window_size = 220 # length of each sequence
step_size = 200 # distance between the start positions of adjacent sequences

# Read in the FASTA file
records = SeqIO.parse(input_file, "fasta")

# Slide the window across each sequence and extract sub-sequences and overhangs
subsequences = []
order = []
for record in records:
    sequence = str(record.seq)
    record_order = []
    
    # Slide the window across the sequence and extract sub-sequences
    prev_subseq_end = None
    for i in range(0, len(sequence)-window_size+1, step_size):
        # Determine the starting position of the current sub-sequence
        if prev_subseq_end is None:
            # First sub-sequence of the sequence
            start_pos = i
        else:
            # Subsequent sub-sequence
            start_pos = prev_subseq_end - 20

        # Extract the current sub-sequence
        subseq = sequence[start_pos:start_pos+window_size]
        subseq_record = SeqRecord(Seq(subseq), id=f"{record.id}_subsequence_{i//step_size+1}", description="")
        subsequences.append(subseq_record)
        record_order.append(len(subsequences)-1)

        # Update the end position of the current sub-sequence
        prev_subseq_end = start_pos + window_size

    overhang_end = sequence[-window_size:]
    if len(overhang_end) == window_size:
        overhang_end_record = SeqRecord(Seq(overhang_end), id=f"{record.id}_overhang_end", description="")
        subsequences.append(overhang_end_record)
        record_order.append(len(subsequences)-1)

    # Add the record order to the global order list
    order.extend(record_order)


# Write the sub-sequences and overhangs to a new FASTA file in the specified order
ordered_subsequences = [subsequences[i] for i in order]
SeqIO.write(subsequences, output_file, "fasta")
