from Bio import Entrez, SeqIO
import time

Entrez.email = "tucorreo@ejemplo.com"  # ← reemplaza con tu correo

accessions = [
    "NZ_CP187663.1",
    "NZ_OX460973.1",
    "NZ_CP128184.1",
    "NZ_CP082243.1",
    "NZ_LN999829.1"
]

for acc in accessions:
    try:
        handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        species = record.annotations.get("organism", "Unknown")
        print(f"{acc}\t{species.replace(' ', '_')}")
        time.sleep(0.4)
    except Exception as e:
        print(f"{acc}\tUnknown")

