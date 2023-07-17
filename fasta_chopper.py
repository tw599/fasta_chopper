import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Extract subsequences and add barcodes.")
parser.add_argument("input_file", help="Input FASTA file")
parser.add_argument("output_file", help="Output FASTA file")
parser.add_argument("barcode_file", help="Barcode file in FASTA format")
parser.add_argument("site1", help="Primer site 1 sequence")
parser.add_argument("site2", help="Primer site 2 sequence")
parser.add_argument(
    "--window_size",
    type=int,
    default=300,
    help="Length of each subsequence (default: 300)",
)
parser.add_argument(
    "--step_size",
    type=int,
    default=275,
    help="Distance between the start positions of adjacent subsequences (default: 275)",
)

args = parser.parse_args()

# Read input FASTA file
records = SeqIO.parse(args.input_file, "fasta")

# Read barcode sequences from the barcode file
barcode_sequences = []
with open(args.barcode_file) as barcode_file:
    for record in SeqIO.parse(barcode_file, "fasta"):
        barcode_sequences.append(str(record.seq))

# Slide the window across each sequence and extract subsequences and overhangs
subsequences = []
order = []
record_counter = 1

for record in records:
    sequence = str(record.seq)
    record_order = []

    # Slide the window across the sequence and extract subsequences
    prev_subseq_end = None
    subseq_counter = 1
    for i in range(0, len(sequence) - args.window_size + 1, args.step_size):
        # Determine the starting position of the current sub-sequence
        if prev_subseq_end is None:
            # First sub-sequence of the sequence
            start_pos = i
        else:
            # Subsequent sub-sequence
            start_pos = prev_subseq_end - 20

        # Extract the current sub-sequence
        subseq = sequence[start_pos:start_pos + args.window_size]

        # Update the ID of the sub-sequence with the serial number
        subseq_id = f"{record.id}_{subseq_counter}"
        subseq_record = SeqRecord(
            Seq(subseq),
            id=subseq_id,
            description="",
        )
        subsequences.append(subseq_record)
        record_order.append(len(subsequences) - 1)

        # Update the sub-sequence counter
        subseq_counter += 1

        # Update the end position of the current sub-sequence
        prev_subseq_end = start_pos + args.window_size

    overhang_end = sequence[-args.window_size:]
    if len(overhang_end) == args.window_size:
        # Update the ID of the overhang_end sequence with the serial number
        overhang_id = f"{record.id}_{subseq_counter}"
        overhang_end_record = SeqRecord(
            Seq(overhang_end),
            id=overhang_id,
            description="",
        )
        subsequences.append(overhang_end_record)
        record_order.append(len(subsequences) - 1)

    # Add the record order to the global order list
    order.extend(record_order)

    # Update the record counter
    record_counter += 1

# Read barcode sequences from the barcode file
barcode_sequences = []
with open(args.barcode_file) as barcode_file:
    for record in SeqIO.parse(barcode_file, "fasta"):
        barcode_sequences.append(str(record.seq))

# Slide the window across each sequence and extract subsequences and overhangs
final_subsequences = []
barcode_index = 0

for subseq_record in subsequences:
    sequence = str(subseq_record.seq)

    # Get the barcode for this sub-sequence
    barcode = barcode_sequences[barcode_index]
    barcode_index = (barcode_index + 1) % len(barcode_sequences)

    # Combine the barcode and sub-sequence
    subseq_with_barcode = barcode + sequence

    # Add primer sites to the 5' and 3' ends
    final_sequence = args.site1 + subseq_with_barcode + args.site2

    subseq_record.seq = Seq(final_sequence)
    final_subsequences.append(subseq_record)

# Write the sub-sequences and overhangs to a new FASTA file in the specified order
ordered_subsequences = [final_subsequences[i] for i in order]
SeqIO.write(ordered_subsequences, args.output_file, "fasta")
