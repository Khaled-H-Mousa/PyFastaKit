
from Module.PySeq import read_fasta , calculate_seq_length, verify_dna_bases, find_start_codon, find_stop_codon, DNA_Transcription, DNA_Translate

FastaFile = "sequence.fasta"

myseqs = read_fasta(FastaFile)
mylenght = calculate_seq_length(myseqs)
NonDNA = verify_dna_bases(myseqs)
STARTCodon = find_start_codon(myseqs)
ENDCodon = find_stop_codon(myseqs)
dnaTrans = DNA_Transcription(myseqs)
dnatranslate = DNA_Translate(myseqs)
print(ENDCodon)
#print(dnatranslate)
#print(STARTCodon)
#print(myseqs)
#print(NonDNA)
#print(mylenght)