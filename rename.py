#!/usr/bin/env python3
import csv
import argparse
import re

def load_mapping(csv_file):
    mapping = {}
    with open(csv_file, newline='') as f:
        sample = f.read(1024)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t'])
        except csv.Error:
            dialect = csv.get_dialect('excel')
        reader = csv.reader(f, dialect)
        for row in reader:
            if len(row) < 2:
                continue
            mapping[row[0].strip()] = row[1].strip()
    return mapping

def rename_fasta(in_fasta, mapping, out_fasta, gene_name):
    with open(in_fasta) as fin, open(out_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                core = header.split('|', 1)[1] if '|' in header else header
                m = re.match(r"(.+?)_cds", core)
                if m:
                    seqid = m.group(1)
                else:
                    seqid = core.split()[0]
                desc = mapping.get(seqid, '')
                parts = [seqid]
                if desc:
                    parts.append(desc)
                if gene_name:
                    parts.append(gene_name)
                fout.write(">" + "\t".join(parts) + "\n")
            else:
                fout.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renombra FASTA usando un CSV de mapeo e inserta nombre de gen")
    parser.add_argument('-i', '--input',  required=True, help='FASTA de entrada')
    parser.add_argument('-m', '--mapping', required=True, help='CSV de mapeo (ID<tab>Descripción)')
    parser.add_argument('-o', '--output', required=True, help='FASTA de salida')
    args = parser.parse_args()

    # Pregunta interactiva por el nombre del gen
    gene_name = input("Introduce el nombre del gen para incluir en el encabezado: ").strip()
    mapping = load_mapping(args.mapping)
    rename_fasta(args.input, mapping, args.output, gene_name)
