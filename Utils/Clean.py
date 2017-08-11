import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Returns a fasta with aligned subsequences', add_help=True)
parser.add_argument('-t','--tabular', dest='tabular', metavar='inputFile', type=str, help='File containing start and end position', required=True)
parser.add_argument('-f','--fasta', dest='genomeFasta', metavar='inputFile', type=str, help='Fasta containing the entire sequence', required=True)
parser.add_argument('-ti','--taxid', dest='taxid', metavar='organism id', type=int, help='organism id in ncbi', required=True)

args = parser.parse_args()

InputFile = open(args.tabular, "r")
OutputFile = open(args.tabular + "-Cured.txt", "w")
LogFile = open(args.tabular + "-log.txt", "w")

Iterator = SeqIO.parse(args.genomeFasta, "fasta")

# Checking for duplicated fasta entries

seqsFasta=[]

for record in Iterator:
	if record.id in seqsFasta:
		print("Duplicated sequence in genome FASTA")
		os.exit()
	else
		seqsFasta.append(record.id)

# Secondaries and discontinued

for line in InputFile:
	splitList = line.split('\t')

	if(not splitList[0] == int(args.taxid)):
		LogFile.write(str(splitList[3]) + ": Gene ID is not from the selected species\n")
		continue

	if(splitList[4] == "discontinued"):
		if (not splitlist[11] in seqsFasta):
			LogFile.write(str(splitList[11]) + ": Nucleotide identifier for gene " +  str(splitList[2]) + " not in input FASTA\n")
		LogFile.write(str(splitList[3]) + ": Discontinued gene identifier\n")
	elif(splitList[4] == "secondary"):
		if (not splitlist[11] in seqsFasta):
			LogFile.write(str(splitList[11]) + ": Nucleotide identifier for gene " +  str(splitList[3]) + " not in input FASTA\n")
		if(str(splitList[12]) != '' and str(splitList[13]) != ''):
			OutputFile.write(line)
		else:
			LogFile.write(str(splitList[3]) + ": No start or end positions [secondary]\n")
	else:
		if(str(splitList[12]) != '' and str(splitList[13]) != ''):
			OutputFile.write(line)
		else:
			LogFile.write(str(splitList[2]) + ": No start or end positions [live]\n")
		if (not splitlist[11] in seqsFasta):
			LogFile.write(str(splitList[11]) + ": Nucleotide identifier for gene " +  str(splitList[2]) + " not in input FASTA\n")


InputFile.close()
OutputFile.close()
LogFile.close()