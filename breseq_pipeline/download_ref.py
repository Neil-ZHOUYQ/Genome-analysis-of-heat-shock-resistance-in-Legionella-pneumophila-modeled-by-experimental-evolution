# This script uses biopython to download reference genome from ncbi database
# Author: Jason
# Need to install biopython first
from Bio import Entrez, SeqIO

# Email address
Entrez.email = "lwy020728@gmail.com"  

# Input NO.
accession_number = str(input("Please enter the ref_number: "))

try:
    # Download FASTA file
    fasta_handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
    fasta_record = fasta_handle.read()
    
    with open(f"{accession_number}.fasta", "w") as fasta_output:
        fasta_output.write(fasta_record)

    print(f"FASTA Download Complete: {accession_number}.fasta")

    # Download GFF3 file
    gff3_handle = Entrez.efetch(db="genome", id=accession_number, rettype="gff3", retmode="text")
    gff3_record = gff3_handle.read()
    
    with open(f"{accession_number}.gff3", "w") as gff3_output:
        gff3_output.write(gff3_record)

    print(f"GFF3 Download Complete: {accession_number}.gff3")

except Exception as e:
    print(f"Error: {e}")
