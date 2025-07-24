#!/usr/bin/python3
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


try:
    fasta_file = sys.argv[1]
    listid = sys.argv[2]
except:
    print('ERROR: ./getSeqById.py fasta-file list')
    sys.exit(0)

sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
filename = os.path.basename(fasta_file)
basename = filename.split('.')
seq_records = []

for line in open( listid ):
    id = line.strip('\n')
    if id in sequences:
        print(f'Obteniendo secuencia {id}')
        seq_record = SeqRecord(Seq(str(sequences[id].seq)), id=sequences[id].id, description='')
        seq_records.append(seq_record)
    else:
        print(f'{id} no encontrado')

# write to file
SeqIO.write(seq_records, f'{basename[0]}_genes.fasta', 'fasta')
