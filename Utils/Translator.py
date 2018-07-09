from Bio import SeqIO
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

Data = SeqIO.parse("mirna_405.fasta", "fasta")
RNA = SeqIO.parse("gene_sequence_405.fasta", "fasta")
i=0
Test = []
TestR = []

for record in RNA:
	TestR.append(record)

for record in Data:
	record.seq = record.seq.back_transcribe()
	Test.append(record)

for record in Test:
	i = i + 1
	TestR.append(record) 	
	SeqIO.write(TestR, "AlignmentInput" + str(i) + ".fasta", "fasta")
	TestR.pop(-1)