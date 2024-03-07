

RNA_CODON_LIST = { 
    "U": {
        "U": {
            "U": "Phe",
            "C": "Phe",
            "A": "Leu",
            "G": "Leu",
        },
        "C": {
            "U": "Ser",
            "C": "Ser",
            "A": "Ser",
            "G": "Ser",
        },
        "A": {
            "U": "Tyr",
            "C": "Tyr",
            "A": "STOP",
            "G": "STOP",
        },
        "G": {
            "U": "Cys",
            "C": "Cys",
            "A": "STOP",
            "G": "Trp",
        },
    }, 
    "C": {
        "U": {
            "U": "Leu",
            "C": "Leu",
            "A": "Leu",
            "G": "Leu",
        },
        "C": {
            "U": "Pro",
            "C": "Pro",
            "A": "Pro",
            "G": "Pro",
        },
        "A": {
            "U": "His",
            "C": "His",
            "A": "Gln",
            "G": "Gln",
        },
        "G": {
            "U": "Arg",
            "C": "Arg",
            "A": "Arg",
            "G": "Arg",
        },
    }, 
    "A": {
        "U": {
            "U": "Ile",
            "C": "Ile",
            "A": "Ile",
            "G": ["START", "Met"],
        },
        "C": {
            "U": "Thr",
            "C": "Thr",
            "A": "Thr",
            "G": "Thr",
        },
        "A": {
            "U": "Asn",
            "C": "Asn",
            "A": "Lys",
            "G": "Lys",
        },
        "G": {
            "U": "Ser",
            "C": "Ser",
            "A": "Arg",
            "G": "Arg",
        },
    }, 
    "G": {
        "U": {
            "U": "Val",
            "C": "Val",
            "A": "Val",
            "G": "Val",
        },
        "C": {
            "U": "Ala",
            "C": "Ala",
            "A": "Ala",
            "G": "Ala",
        },
        "A": {
            "U": "Asp",
            "C": "Asp",
            "A": "Glu",
            "G": "Glu",
        },
        "G": {
            "U": "Gly",
            "C": "Gly",
            "A": "Gly",
            "G": "Gly",
        },
    }, 
}

NUCLEOBASES = {
    
    "A": "Adenin",
    "T": "Thymin",
    "C": "Cytosin",
    "G": "Guanin",
    
    "U": "Uracil"
    
}

NUCLEOBASE_PAIR_TRANSLATION = {
    
    "A": "U",
    "T": "A",
    
    "C": "G",
    "G": "C"
    
}


AMINO_ACID_TRANSLATIONS = {
    
    "Ala": "Alanin",
    "Arg": "Arginin",
    "Asn": "Asparagin",
    "Asp": "Asperaginsäure",
    "Cys": "Cystein",
    "Gln": "Glutamin",
    "Glu": "Glutaminsäure",
    "Gly": "Glycin",
    "His": "Histidin",
    "Ile": "Isoleucin",
    "Leu": "Leucin",
    "Lys": "Lysin",
    "Met": "Methionin",
    "Phe": "Phenylalanin",
    "Pro": "Prolin",
    "Ser": "Serin",
    "Thr": "Threonin",
    "Trp": "Tryptophan",
    "Tyr": "Tyrosin",
    "Val": "Valin",
    
}


def print_DNA_mRNA_aminoacid_table(DNA_codons: list[str], mRNA_codons: list[str], amino_acids: list[str]):
    
    cellPadding = 1
    
    def getMaxLen(cells: list):
        return max(len(entry) for entry in cells)
    
    def padding():
        return " "*cellPadding
    def fillWhitespace(curLen: int, allItems: list[str]):
        return " " * (getMaxLen(allItems) - curLen)
    
    
    
    def getThreeRows(codonDNA: list[str], codon_mRNA: list[str], amino_acid: str):
        row1 =  "|" + padding() +  codonDNA[0] +             fillWhitespace(len(codonDNA[0]), NUCLEOBASES.keys()) +                    padding() # Kürzel
        row1 += "|" + padding() + NUCLEOBASES[codonDNA[0]] + fillWhitespace(len(NUCLEOBASES[codonDNA[0]]), NUCLEOBASES.values()) +     padding() # Name
        row1 += "|" + (padding() * 2) + (" " * 2) + (padding() * 2)                                                                              # Free arrow space
        row1 += "|" + padding() + codon_mRNA[0] + fillWhitespace(len(codon_mRNA[0]), NUCLEOBASES.keys()) +                             padding() # Kürzel
        row1 += "|" + padding() + NUCLEOBASES[codon_mRNA[0]] + fillWhitespace(len(NUCLEOBASES[codon_mRNA[0]]), NUCLEOBASES.values()) + padding() # Name
        row1 += "|" + (padding() * 2) + " " + (padding() * 2)                                                                                    # Free '=' space
        row1 += "|" + padding() + fillWhitespace(0, AMINO_ACID_TRANSLATIONS.keys()) + padding()                                                  # Free amino acid kürzel space
        row1 += " " + padding() + fillWhitespace(0, AMINO_ACID_TRANSLATIONS.values()) + padding()                                                # amino acid name free space
        row1 += "|"
        
        row2 =  "|" + padding() +  codonDNA[1] +             fillWhitespace(len(codonDNA[1]), NUCLEOBASES.keys()) +                    padding() # Kürzel
        row2 += "|" + padding() + NUCLEOBASES[codonDNA[1]] + fillWhitespace(len(NUCLEOBASES[codonDNA[1]]), NUCLEOBASES.values()) +     padding() # Name
        row2 += "|" + (padding() * 2) + "=>" + (padding() * 2)                                                                                   # '=>' arrow
        row2 += "|" + padding() + codon_mRNA[1] + fillWhitespace(len(codon_mRNA[1]), NUCLEOBASES.keys()) +                             padding() # Kürzel
        row2 += "|" + padding() + NUCLEOBASES[codon_mRNA[1]] + fillWhitespace(len(NUCLEOBASES[codon_mRNA[1]]), NUCLEOBASES.values()) + padding() # Name
        row2 += "|" + (padding() * 2) + "=" + (padding() * 2)                                                                                    # '='
        row2 += "|" + padding() + amino_acid + fillWhitespace(len(amino_acid), AMINO_ACID_TRANSLATIONS.keys()) + padding()                           # Free amino acid kürzel space
        row2 += " " + padding() + AMINO_ACID_TRANSLATIONS[amino_acid] + fillWhitespace(len(AMINO_ACID_TRANSLATIONS[amino_acid]), AMINO_ACID_TRANSLATIONS.values()) + padding()  # amino acid name
        row2 += "|"
        
        row3 =  "|" + padding() +  codonDNA[2] +             fillWhitespace(len(codonDNA[2]), NUCLEOBASES.keys()) +                    padding() # Kürzel
        row3 += "|" + padding() + NUCLEOBASES[codonDNA[2]] + fillWhitespace(len(NUCLEOBASES[codonDNA[2]]), NUCLEOBASES.values()) +     padding() # Name
        row3 += "|" + (padding() * 2) + (" " * 2) + (padding() * 2)                                                                              # Free arrow space
        row3 += "|" + padding() + codon_mRNA[2] + fillWhitespace(len(codon_mRNA[2]), NUCLEOBASES.keys()) +                             padding() # Kürzel
        row3 += "|" + padding() + NUCLEOBASES[codon_mRNA[2]] + fillWhitespace(len(NUCLEOBASES[codon_mRNA[2]]), NUCLEOBASES.values()) + padding() # Name
        row3 += "|" + (padding() * 2) + " " + (padding() * 2)                                                                                    # Free '=' space
        row3 += "|" + padding() + fillWhitespace(0, AMINO_ACID_TRANSLATIONS.keys()) + padding()                                                  # Free amino acid kürzel space
        row3 += " " + padding() + fillWhitespace(0, AMINO_ACID_TRANSLATIONS.values()) + padding()                                                # amino acid name free space
        row3 += "|"
        
        return [row1, row2, row3]
    
    
    def getAllRows():
        out = []
        allRows = []
        for i in range(len(DNA_codons)):
            rows = getThreeRows(DNA_codons[i], mRNA_codons[i], amino_acids[i])
            allRows.append(rows)
        
        for i in range(len(allRows)):
            
            out.append("-" * len(allRows[0][0]))
            out.append(allRows[i][0])
            out.append(allRows[i][1])
            out.append(allRows[i][2])
        
        out.append("-" * len(allRows[0][0]))
        return out
    
    
    
    for line in getAllRows():
        print(line)
            
    
    



def translateDNA_to_mRNA(DNA: list[str]):
    return [NUCLEOBASE_PAIR_TRANSLATION[base] for base in DNA]


def process_DNA(DNA: list[str]):
    dna_codon_list: list[str] = []
    mrna_codon_list: list[str] = []
    amino_acids: list[str] = []
    
    codon_reader: list[str] = [] # max len 3
    is_reading: bool = False
    
    for base in DNA:
        if len(codon_reader) == 3:
            mrnaCodonList = translateDNA_to_mRNA(codon_reader)
            codonDNA = "".join(codon_reader)
            codonMRNA = "".join(mrnaCodonList)
            amino_acid = RNA_CODON_LIST[mrnaCodonList[0]][mrnaCodonList[1]][mrnaCodonList[2]]
            
            if isinstance(amino_acid, list):
                if "START" in amino_acid and not is_reading:
                    is_reading = True
                    codon_reader = [base]
                    continue
                elif "Met" in amino_acid and is_reading:
                    amino_acid = "Met"
                    
            if amino_acid == "STOP":
                return [dna_codon_list, mrna_codon_list, amino_acids]
            
            if is_reading:
                amino_acids.append(amino_acid)
                dna_codon_list.append(codonDNA)
                mrna_codon_list.append(codonMRNA)
                
            codon_reader = []
            
        codon_reader.append(base)
            
    return [dna_codon_list, mrna_codon_list, amino_acids]
  
        
def get_amino_acids_from_dna(DNA: list[str]):
    result = process_DNA(DNA)
    
    
    print_DNA_mRNA_aminoacid_table(result[0], result[1], result[2])
    
    
get_amino_acids_from_dna(list("AACTACGATGACTATGGAGAACCCACTCGGGTTGCACTAAA"))


