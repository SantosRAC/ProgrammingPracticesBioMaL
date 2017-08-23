# 
__author__      = "Henrique Frajacomo"

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

parser = argparse.ArgumentParser(description='Returns a fasta with aligned subsequences', add_help=True)
parser.add_argument('-t','--tabular', dest='tabular', metavar='inputFile', type=str, help='File containing start and end position', required=True)
parser.add_argument('-f','--fasta', dest='genomeFasta', metavar='inputFile', type=str, help='Fasta containing the entire sequence', required=True)
parser.add_argument('-o','--out', dest='geneFasta', metavar='outputFile', type=str, help='Fasta containing gene sequences', required=True)

args = parser.parse_args()
argTable = args.tabular
sequenceFile = args.genomeFasta
resultFile = args.geneFasta
geneIdCount = {}

resultFileOBJ = open(resultFile, "w")

Iterator = SeqIO.parse(sequenceFile, "fasta")
numberOfGenes=0
numberOfSeqs=0
for record in Iterator:
	numberOfSeqs=numberOfSeqs+1
	dataTable = open(argTable)
	for gene in dataTable:
		splitList = gene.split('\t')
		if(splitList[11] == record.id):
			start = int(splitList[12])-1
			end = int(splitList[13])
			numberOfGenes=numberOfGenes+1
			seq_description="[start:" + splitList[12] + "] [end:" + splitList[13] + "] [chrom:" + splitList[10] + "] [strand:" + splitList[14] +"]"
			if(splitList[14]=="plus"):
				generecord = SeqRecord(Seq(str(record.seq[start:end])),id=splitList[2],description=seq_description)
				if generecord.id in geneIdCount.keys():
					geneIdCount[str(generecord.id)]=geneIdCount[str(generecord.id)]+1
				else:
					geneIdCount[str(generecord.id)]=1
				generecord.id=str(generecord.id) + "_" + str(geneIdCount[str(generecord.id)])
				SeqIO.write(generecord,resultFileOBJ,"fasta")
			else:
				sequence=Seq(str(record.seq[start:end]))
				generecord = SeqRecord(sequence.reverse_complement(),id=splitList[2],description=seq_description)
				if generecord.id in geneIdCount.keys():
					geneIdCount[str(generecord.id)]=geneIdCount[str(generecord.id)]+1
				else:
					geneIdCount[str(generecord.id)]=1
				generecord.id=str(generecord.id) + "_" + str(geneIdCount[str(generecord.id)])
				SeqIO.write(generecord,resultFileOBJ,"fasta")
	dataTable.close()

print("There are " + str(numberOfGenes) + " genes in the input TAB file! Yeahhh")
print("There are " + str(numberOfSeqs) + " sequences in the FASTA file! Yeahhh")

Iterator.close()
resultFileOBJ.close()
