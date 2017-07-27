import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Returns a fasta with aligned subsequences', add_help=True)
parser.add_argument('-t','--tabular', dest='tabular', metavar='inputFile', type=str, help='File containing start and end position', required=True)
parser.add_argument('-f','--fasta', dest='fasta', metavar='inputFile', type=str, help='Fasta containing the entire sequence', required=True)

args = parser.parse_args()
argTable = args.tabular
sequenceFile = args.fasta


repeat = False
fullList = []
splitList = []
dataTable = open(argTable)
resultFile = open("ResultFile.fasta", "w")

for line in dataTable:
	if(repeat):
 		fullList.append(line.strip('\n'))
 	repeat = True

repeat = False

StartIterator = SeqIO.parse(sequenceFile, "fasta")
Iterator = SeqIO.parse(sequenceFile, "fasta")

for element in fullList:
	Iterator = StartIterator
	splitList = element.split('\t')
	for record in Iterator:
		if(splitList[11] == record.id):
			start = int(splitList[12])
			end = int(splitList[13])
			resultFile.write(">" + splitList[2] +"\n" + str(record.seq[start:end]))
			print(end-start)


dataTable.close()
resultFile.close()

'''
for record in SeqIO.parse("sequence.fasta", "fasta"):
	if record.id in fullList:
		resultFile.write(">" + record.id)
		start = 
'''



 #print(">%s\n%s\n%d" % (record.id, record.seq, len(record)))
