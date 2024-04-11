# iupac_codes.py
from Bio.Data import CodonTable


# NUCLEOTIDE DICTIONARIES
Standard_Nucleotides = {
    "A": "adenosine",
    "C": "cytidine",
    "G": "guanine",
    "T": "thymidine",
    "N": "any (A/G/C/T)",
    "U": "uridine",
}

Degenerate_Nucleotides = {
    "K": "keto (G/T)",
    "S": "strong (G/C)",
    "Y": "pyrimidine (T/C)",
    "M": "amino (A/C)",
    "W": "weak (A/T)",
    "R": "purine (G/A)",
    "B": "G/T/C",
    "D": "G/A/T",
    "H": "A/C/T",
    "V": "G/C/A",
    "-": "gap of indeterminate length",
}

Nucleotide_Complements = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "U": "A",
    "K": "M",
    "S": "S",
    "Y": "R",
    "M": "K",
    "W": "W",
    "R": "Y",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "-": "-",
}

# AMINOACID DICTIONARIES
Standard_AminoAcids = {
    "A": "alanine",
    "B": "aspartate/asparagine",
    "C": "cystine",
    "D": "aspartate",
    "E": "glutamate",
    "F": "phenylalanine",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "K": "lysine",
    "L": "leucine",
    "M": "methionine",
    "N": "asparagine",
    "P": "proline",
    "Q": "glutamine",
    "R": "arginine",
    "S": "serine",
    "T": "threonine",
    "U": "selenocysteine",
    "V": "valine",
    "W": "tryptophan",
    "Y": "tyrosine",
    "Z": "glutamate/glutamine",
    "X": "any",
    "*": "translation stop",
}

# CODON TABLES
# Access the standard DNA codon table (ID 1)
standard_table = CodonTable.unambiguous_dna_by_id[1]
