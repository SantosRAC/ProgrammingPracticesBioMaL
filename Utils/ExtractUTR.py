import Bio
import numpy
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Extracts CDS and mRNA entries from a gff3', add_help=True)
parser.add_argument('-g','--gff', dest='gff', metavar='inputGFF', type=str, help='GFF containing information about features in genome', required=True)

args = parser.parse_args()

start_time = datetime.now()

'''
file = open(args.gff, "r")

mrnafile = open("mrna.gff", "w")
cdsfile = open("cds.gff", "w")

print("Creating CDS and mRNA separate files...")

# Creates a file containing only mRNA and CDS entries ---
for line in file:
	splitList = line.split('\t')
	if(line[0] == "#"):
		continue
	elif(splitList[2] == "mRNA"):
		mrnafile.write(line)
	elif(splitList[2] == "CDS"):
		cdsfile.write(line)

file.close()
mrnafile.close()
cdsfile.close()
# -------------------------------------------------------
'''

# Sorts mRNA file
'''
print("Starting mRNA file sorting...")
mrnafile = open("mrna.gff", "r")
mrna = []
numbers = []
linha = mrnafile.readline()

while linha:  # Fills all geneID Entries into a list
	mrna.append(linha)
	splitmRNA = linha.split("\t")
	splitmRNAData = splitmRNA[8].split(";")

	CompareID = (((splitmRNAData[2].split("="))[1].split(","))[0].split(":"))[1]
	numbers.append(int(CompareID))

	linha = mrnafile.readline()

mrnafile.close()
numbers.sort()

print("Sorting memory elements...")
finalList = []

last_num = 0
last_pos = 0
i = 0

while numbers:
	i = last_pos
	while(i < len(mrna)):
		splitmRNA = mrna[i].split("\t")
		splitmRNAData = splitmRNA[8].split(";")
		CompareID = (((splitmRNAData[2].split("="))[1].split(","))[0].split(":"))[1]
		if(int(CompareID) == numbers[0]):
			finalList.append(mrna[i])
			mrna.pop(i)
			try:
				if(numbers[1] == numbers[0]):  # Partitioning
					last_pos = i
			except IndexError:
				pass
			else:
				last_pos = 0
			last_num = numbers[0]
			numbers.pop(0)
			i = i + 1
			break
		i = i + 1

	if(len(numbers)%100 == 0):
		print(len(numbers))

print("Time Elapsed: " + str(datetime.now() - start_time))

mrnafile = open("mrna.gff", "w")
print("Writing sorted mRNA file...")

for element in finalList:
	mrnafile.write(element)

logfile = open("timelog.txt", "w")
print("Time Elapsed: " + str(datetime.now() - start_time))
logfile.write("Time Taken: " + str(datetime.now() - start_time))

mrnafile.close()
finalList = []
'''
# ------------------------------------------------------
# Sorts CDS file

print("Starting CDS file sorting...")
cdsfile = open("cds.gff", "r")
cds = []
numbers = []
linha = cdsfile.readline()

while linha:  # Fills all geneID Entries into a list
	cds.append(linha)
	splitCDS = linha.split("\t")
	splitCDSData = splitCDS[8].split(";")

	if(((splitCDSData[2].split("="))[1].split(":"))[0] == "CCDS" or ((splitCDSData[2].split("="))[1].split(":"))[0] == "Genbank"):
		CompareID = (((splitCDSData[2].split("="))[1].split(","))[1].split(":"))[1]
	else:
		CompareID = (((splitCDSData[2].split("="))[1].split(","))[0].split(":"))[1]

	numbers.append(int(CompareID))

	linha = cdsfile.readline()

cdsfile.close()
numbers.sort()

print("Sorting memory elements...")
finalList = []

last_num = 0
last_pos = 0
i = 0

while numbers:
	i = last_pos
	while(i < len(cds)):
		splitCDS = cds[i].split("\t")
		splitCDSData = splitCDS[8].split(";")

		if(((splitCDSData[2].split("="))[1].split(":"))[0] == "CCDS" or ((splitCDSData[2].split("="))[1].split(":"))[0] == "Genbank"):
			CompareID = (((splitCDSData[2].split("="))[1].split(","))[1].split(":"))[1]
		else:
			CompareID = (((splitCDSData[2].split("="))[1].split(","))[0].split(":"))[1]

		if(int(CompareID) == numbers[0]):
			finalList.append(cds[i])
			cds.pop(i)
			try:
				if(numbers[1] == numbers[0]):  # Partitioning
					last_pos = i
			except IndexError:
				pass
			else:
				last_pos = 0
			last_num = numbers[0]
			numbers.pop(0)
			i = i + 1
			break
		i = i + 1

	if(len(numbers)%100 == 0):
		print(len(numbers))

print("Time Elapsed: " + str(datetime.now() - start_time))


mrnafile = open("cds.gff", "w")
print("Writing sorted CDS file...")

for element in finalList:
	mrnafile.write(element)

logfile = open("timelog.txt", "a")
logfile.write("\n")
print("Time Elapsed: " + str(datetime.now() - start_time))
logfile.write("Time Taken: " + str(datetime.now() - start_time))

mrnafile.close()
finalList = []
'''
#-----------------------------------------------------------
'''
'''
# ------------- Unused ----------------------------------------------
mrnafile = open("mrna.gff", "r")

logfile = open("cleaninglog.txt", "w")
logfile.close()

linecount = 0


# Creates an line offset list for cds iteration
with open("cds.gff", "r") as f:
	cdsline = []
	offset = 0
	linha = f.readline()

	while linha:
		offset = f.tell()
		cdsline.append(offset)
		linha = f.readline()

'''
# ---------------------------------------------------------------
'''	
    

# Create partitions on the CDS file for quick indexing ----
print("Creating an entry index...")
lastelement = ''

IDlist = [0]
with open("cds.gff", "r") as f:
	linha = f.readline()

	while linha:
		if((((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS"):
			line_element = ((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1]
		else:
			line_element = ((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0]
		if(lastelement == ''):
			lastelement = line_element
		elif(line_element == lastelement):
			linha = f.readline()
			continue
		else:
			IDlist.append(f.tell())
			lastelement = line_element
		linha = f.readline()
print("Size: " + str(len(IDlist)))

print("Time Elapsed: " + str(datetime.now() - start_time))

#------------------------------------------------------------
# Starts ORF extraction

cdsfile = open("cds.gff", "r")
lastGeneID = ''
IDIndex = 0
resultfile = open("mrnainfo.gff", "w")
resultfile.write("#GeneID" + "\t" + "rna_code" + "\t" + "lowest_entry" + "\t" + "highest_entry" + "\n")

resultfile.close()

print("Starting mRNAs' ORF extraction...")

for line in mrnafile:
	highestEntry = 0
	lowestEntry = 900000000000000000000
	print("Processing mRNA #" + str(linecount + 1))
	linecount = linecount + 1

	splitmRNA = line.split("\t")
	splitmRNAData = splitmRNA[8].split(";")
	if(lastGeneID == ''):
		lastGeneID = ((splitmRNAData[2].split("="))[1].split(","))[0]
	elif(lastGeneID != ((splitmRNAData[2].split("="))[1].split(","))[0]):
		lastGeneID = ((splitmRNAData[2].split("="))[1].split(","))[0]
		IDIndex = IDIndex + 1
		cdsfile.seek(IDlist[IDIndex],0)
	else:
		cdsfile.seek(IDlist[IDIndex],0)

	found = False
	for line in cdsfile:
		splitCDS = line.split("\t")
		splitCDSData = splitCDS[8].split(";")
		if(((splitCDSData[2].split("="))[1].split(":"))[0] == "CCDS"):
			CompareID = ((splitCDSData[2].split("="))[1].split(","))[1]
		else:
			CompareID = ((splitCDSData[2].split("="))[1].split(","))[0]
		if(CompareID != lastGeneID):  # stop iterating if out of index range
			print(str(lastGeneID) + "\t" + str(CompareID))
			break
		if(lastGeneID == CompareID and (splitmRNAData[0].split("="))[1] == (splitCDSData[1].split("="))[1]):
			found = True
			if(int(splitCDS[3]) < lowestEntry):
				lowestEntry = int(splitCDS[3])
			if(int(splitCDS[4]) > highestEntry):
				highestEntry = int(splitCDS[4])
	if(not found):
		logfile = open("cleaninglog.txt", "a")
		logfile.write(((splitmRNAData[2].split("="))[1].split(","))[0] + "\t" + (splitmRNAData[0].split("="))[1] + "\t" + "CDS not found\n")
		logfile.close()
	else:
		resultfile = open("mrnainfo.gff", "a")
		resultfile.write(((splitmRNAData[2].split("="))[1].split(","))[0] + "\t" + (splitmRNAData[0].split("="))[1] + "\t" + str(lowestEntry) + "\t" + str(highestEntry) + "\n")
		resultfile.close()

print("Process Finished!")
print("Total Time Elapsed: " + str(datetime.now() - start_time))
'''
