
#FastaFile = "/home/khaled/Documents/Python_Tasks/sequence.fasta"
import os
import sys
import re

def read_fasta(FastaFile):
    seqs = {}
    with open(FastaFile) as f:
        header = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
                newheader = header.replace(">", "")
                #header = line[1:]  # Remove ">" from header
                newheader = " ".join(newheader.split(" ")[:3])
                seqs[newheader] = ""
            else:
                seqs[newheader] += line.upper()
                #seqs[newheader] += line.lower()
    return seqs



def calculate_seq_length(seqs):
    """
    Calculates the lengths of multiple sequences.
    """         
    get_length = {}
    for seq_id, seq in seqs.items():
        # calculate sequence length
        get_length[seq_id] = len(seq)
    return get_length



def verify_dna_bases(seqs):
    """
    Verifies that DNA sequences contain only valid nucleotide bases (A, C, G, T).
    """
    for seq_id, seq in seqs.items():
        for nucleotide in seq:
            # check if DNA sequence had non-DNA nucleotide
            if nucleotide not in "ACGT":
                print(f"Error: {seq_id} contains non-DNA nucleotide {nucleotide}")
                return True  # if you want to exit on the first error found
    return False
    #return True



def find_start_codon(seqs):
    """
    Finds the positions of start codons in multiple DNA sequences.
    """
    start_codon = {}
    for seq_id, seq in seqs.items():
        StartCodon = "ATG"
        # Find start codon "ATG"
        start_codon[seq_id] = seq.find(StartCodon)
        #for StartCodon in seq:
            #if StartCodon in seq:
                #print(f"Yes: {header} contain a start codon (ATG)")
    return start_codon
    


#def find_stop_codon(seqs):
#    stop_codon = {}
#    for seq_id, seq in seqs.items():
#        StopCodon = ["TAA", "TAG", "TGA"]
#        NewStopCodon = " ".join(StopCodon)
#        # Find stop codon ("TAA", "TAG", "TGA")
#        stop_codon[seq_id] = seq.find(NewStopCodon)
#    return stop_codon


def find_stop_codon(seqs):
    """
    Finds the positions of stop codons in multiple DNA sequences.
    """
    stop_codon = {}
    for seq_id, seq in seqs.items():
        StopCodon = ["TAA", "TAG", "TGA"]
        Position = []
        for i in range(0,len(seq),3):
            codon = seq[i:i+3]
            if codon in StopCodon:
                Position.append((i,codon))
        stop_codon[seq_id] = Position
    return stop_codon



def DNA_Transcription(seqs):
    """
    Transcribes DNA sequences into RNA by replacing thymine (T) with uracil (U).
    """
    transcripte_dna = {}
    for seq_id, seq in seqs.items():
        # replace "T" base to "U" base
        transcripte_dna[seq_id] = seq.replace("T","U")
    return transcripte_dna



def DNA_Translate(seqs):
    """
    Translates DNA sequences into protein sequences based on codons.
    """
    translate_dna = {}
    # Define the codon to amino acid translation table
    table = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
              'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
              'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
              'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
              'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
              'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
              'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
              'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
              'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
              'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
              'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
              'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
              'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
              'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
              'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
              'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W' }

    # Iterate over each sequence
    #Loops through each DNA sequence in seqs: 
    # seqs.items() returns both the header,sequence1
    for seq_id, seq in seqs.items():
        # Initializes an empty string "protein"
        protein = ""  
        # Translate each codon in the DNA sequence
        # 3 Steps to get codons
        for i in range(0, len(seq), 3):  
            # seq[i:i+3] selects a substring starting at position `i` and ending just before `i+3`
            codon = seq[i:i+3]  # Extract codon
            #The get() method returns the value of the item with the specified key.
            # Translate codon to amino acid, or 'None' if not found
            amino_acid = table.get(codon, "None") 
            protein += amino_acid  # Add amino acid to protein sequence
        
        translate_dna[seq_id] = protein  # Store translated protein for this sequence
    
    return translate_dna

