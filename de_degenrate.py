from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product
import sys

input_file = sys.argv[1]
output_file = sys.argv[1].split('.')[0] + '.degenerate_primer.fas'

degenerate_bases = {
         "R": "AG",
         "Y": "CT",
         "S": "GC",
         "W": "AT",
         "K": "GT",
         "M": "AC",
         "B": "GCT",
         "D": "GAT",
         "H": "ACT",
         "V": "GCA",
         "N": "AGCT"
}

for record in SeqIO.parse(input_file, "fasta"):
    non_degenerate_primers = []
    degenerate_combinations = product(*[degenerate_bases.get(base, base) for base in record.seq])
    c = 1
    for combination in degenerate_combinations:
        non_degenerate_seq = ''.join(combination)
        non_degenerate_primer = SeqRecord(Seq(non_degenerate_seq), id=f"{record.id}_{c}")
        non_degenerate_primers.append(non_degenerate_primer)
        c+=1

    with open(output_file, "w") as f:
        for i in non_degenerate_primers:
             f.write(f">{i.id}\n{i.seq}\n")
