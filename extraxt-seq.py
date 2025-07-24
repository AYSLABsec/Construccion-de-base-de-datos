from Bio import SeqIO
import sys

def extract_sequences(fasta_file, id_file, output_file):
    # Leer los IDs de interés desde el archivo de texto
    with open(id_file, 'r') as id_f:
        ids_to_extract = set(line.strip() for line in id_f)

    # Leer el archivo FASTA y filtrar las secuencias
    with open(fasta_file, 'r') as fasta_f, open(output_file, 'w') as output_f:
        for record in SeqIO.parse(fasta_f, "fasta"):
            if record.id in ids_to_extract:
                SeqIO.write(record, output_f, "fasta")

# Configuración de archivos
fasta_file = sys.argv[1]  # Archivo FASTA de entrada
id_file = sys.argv[2]            # Archivo de texto con los IDs
output_file = "extracted_sequences.fasta"  # Archivo FASTA de salida

# Llamar a la función
extract_sequences(fasta_file, id_file, output_file)
print(f"Secuencias extraídas guardadas en: {output_file}")

