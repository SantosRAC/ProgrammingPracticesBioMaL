# 
__author__      = "Henrique Frajacomo"

import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(description='Returns a fasta with aligned subsequences', add_help=True)
parser.add_argument('-t','--tabular', dest='tabular', metavar='inputFile', type=str, help='File containing start and end position', required=True)
parser.add_argument('-f','--fasta', dest='fasta', metavar='inputFile', type=str, help='Fasta containing the entire sequence', required=True)

args = parser.parse_args()
argTable = args.tabular
sequenceFile = args.fasta

repeat = False
splitList = []

resultFile = open("ResultFile.fasta", "w")

Iterator = SeqIO.parse(sequenceFile, "fasta")
numberOfGenes=0
numberOfSeqs=0
for record in Iterator:
	numberOfSeqs=numberOfSeqs+1
	dataTable = open(argTable)
	for line in dataTable:
		splitList = line.split('\t')
		if(splitList[11] == record.id):
			start = int(splitList[12])
			end = int(splitList[13])
			numberOfGenes=numberOfGenes+1
			if(splitList[14]=="plus"):
				resultFile.write(">" + splitList[2] + "\t" + splitList[14] +"\n" + str(record.seq[start:end]) + "\n")
			else:
				resultFile.write(">" + splitList[2] + "\t" + splitList[14] +"\n" + str((record.seq[start:end])[::-1]) + "\n")  #[::-1] -> reads the entry backwards
	dataTable.close()

print("There are ", numberOfGenes, " genes in the input TAB file! Yeahhh")
print("There are ", numberOfSeqs, " sequences in the FASTA file! Yeahhh")

Iterator.close()
resultFile.close()
