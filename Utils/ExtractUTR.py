import Bio
import numpy
import argparse
from datetime import datetime
import sys
import os
import math

def Partition_mRNA(A, start, end):
	Pivot = int((((((A[end].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1])
	Pindex = start
	i = 0
	for i in range(start,end):
		if(int((((((A[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1]) <= Pivot):
			A[i], A[Pindex] = A[Pindex], A[i]
			Pindex = Pindex + 1
	A[end], A[Pindex] = A[Pindex], A[end]
	return Pindex, A

def QuickSort_mRNA(A, start, end):
	if(start < end):
		Pindex, A = Partition_mRNA(A, start, end)
		A = QuickSort_mRNA(A, start, Pindex-1)
		A = QuickSort_mRNA(A, Pindex+1, end)
	return A

def Partition_CDS(A, start, end):
	if((((((A[end].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((A[end].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
		Pivot = int((((((A[end].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1].split(":"))[1])
	else:
		Pivot = int((((((A[end].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1])
	Pindex = start
	i = 0
	for i in range(start,end):
		if((((((A[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((A[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
			CompareID = (((((A[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1].split(":"))[1]
		else:
			CompareID = (((((A[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1]
		
		if(int(CompareID) <= Pivot):
			A[i], A[Pindex] = A[Pindex], A[i]
			Pindex = Pindex + 1
	A[end], A[Pindex] = A[Pindex], A[end]
	return Pindex, A

def QuickSort_CDS(A, start, end):
	if(start < end):
		Pindex, A = Partition_CDS(A, start, end)
		A = QuickSort_CDS(A, start, Pindex-1)
		A = QuickSort_CDS(A, Pindex+1, end)
	return A


sys.setrecursionlimit(1000000) # Default = 10000
parser = argparse.ArgumentParser(description='Extracts UTR length from mRNA entries', add_help=True)
parser.add_argument('-g','--gff', dest='gff', metavar='inputGFF', type=str, help='GFF containing information about features in genome', required=True)

args = parser.parse_args()

start_time = datetime.now()


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

# Sorts mRNA file

print("Starting mRNA file sorting...")
mrnafile = open("mrna.gff", "r")
mrna = []
linha = mrnafile.readline()

while linha: 
	mrna.append(linha)
	linha = mrnafile.readline()

mrnafile.close()

mrna = QuickSort_mRNA(mrna, 0, len(mrna)-1)


print("Writing sorted file...")
mrnafile = open("mrna.gff", "w")
for element in mrna:
	mrnafile.write(element)
mrnafile.close()


logfile = open("timelog.txt", "w")
print("Time Elapsed: " + str(datetime.now() - start_time))
logfile.write("Time Taken: " + str(datetime.now() - start_time))


# ------------------------------------------------------
# Sorts CDS file

print("Starting CDS file sorting...")
cdsfile = open("cds.gff", "r")
cds = []
linha = cdsfile.readline()
dimension_amount = 0
loops = 0
cdsaux = [[]]
Segmented = False
FileSize = 150000

while linha: 
	cds.append(linha)
	linha = cdsfile.readline()

cdsfile.close()
print("File Size (in lines): " + str(len(cds)))

# File Segmentation if neccessary
if(len(cds)>FileSize):
	Segmented = True
	print("File too big, segmentating it in " + str(math.ceil((len(cds)/(FileSize-1)))) + " partitions")
	while(len(cds)>FileSize):
		for i in range(0,FileSize-1):
			cdsaux[dimension_amount].append(cds.pop(0))
			if(i%10000 == 0):
				print("Analyzed " + str(i) + " lines")
		loops = loops + 1
		dimension_amount = dimension_amount + 1
		cdsaux.append([])
	for i in range(0, len(cds)-1):
		if(i%10000 == 0):
			print("Analyzed " + str(i) + " lines")
		cdsaux[dimension_amount].append(cds.pop(0))

if(Segmented):
	print("Started Sorting Process")
	for i in range(0, dimension_amount-1):
		cdsaux[i] = QuickSort_CDS(cdsaux[i], 0, len(cdsaux[i])-1)

	for i in range(0, dimension_amount-1):
		while(len(cdsaux[i])>0):
			cds.append(cdsaux[i].pop(0))
	cds = QuickSort_CDS(cds, 0, len(cds)-1)

else:
	cds = QuickSort_CDS(cds, 0, len(cds)-1)

print("Writing sorted file...")
cdsfile = open("cds.gff", "w")
for element in cds:
	cdsfile.write(element)
cdsfile.close()


logfile = open("timelog.txt", "a")
print("Time Elapsed: " + str(datetime.now() - start_time))
logfile.write("Time Taken: " + str(datetime.now() - start_time))

# ---------------------------------------------------------------

# Create partitions on the CDS file for quick indexing ----
print("Creating an entry index...")
lastelement = ''

IDlist = [0]
with open("cds.gff", "r") as f:
	linha = f.readline()

	while linha:
		if((((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
			line_element = ((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1]
		else:
			line_element = ((((linha.split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0]
		if(lastelement == ''):
			lastelement = line_element
		elif(line_element != lastelement):
			IDlist.append(last_address)
			lastelement = line_element
			last_address = f.tell()
			continue

		last_address = f.tell()
		linha = f.readline()
print("Size: " + str(len(IDlist)))

print("Time Elapsed: " + str(datetime.now() - start_time))

#------------------------------------------------------------
# Starts ORF extraction

mrnafile = open("mrna.gff", "r")
cdsfile = open("cds.gff", "r")
resultfile = open("mrnainfo.gff", "w")
resultfile.write("#GeneID" + "\t" + "rna_code" + "\t" + "lowest_entry" + "\t" + "highest_entry" + "\n")
resultfile.close()
logfile = open("cleaninglog.txt", "w")
logfile.close()

lastGeneID = ''
linecount = 0
IDIndex = 0

print("Starting mRNAs' ORF extraction...")

linha_2 = mrnafile.readline()

while linha_2:
	highestEntry = 0
	lowestEntry = 900000000000000000000
	print("Processing mRNA #" + str(linecount + 1))
	linecount = linecount + 1

	splitmRNA = linha_2.split("\t")
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
	linha = cdsfile.readline()

	while linha:
		splitCDS = linha.split("\t")
		splitCDSData = splitCDS[8].split(";")
		if(((splitCDSData[2].split("="))[1].split(":"))[0] == "CCDS" or ((splitCDSData[2].split("="))[1].split(":"))[0] == "Genbank"):
			CompareID = ((splitCDSData[2].split("="))[1].split(","))[1]
		else:
			CompareID = ((splitCDSData[2].split("="))[1].split(","))[0]
		if(CompareID != lastGeneID):  # stop iterating if out of index range
			break
		if(lastGeneID == CompareID and (splitmRNAData[0].split("="))[1] == (splitCDSData[1].split("="))[1]):
			found = True
			if(int(splitCDS[3]) < lowestEntry):
				lowestEntry = int(splitCDS[3])
			if(int(splitCDS[4]) > highestEntry):
				highestEntry = int(splitCDS[4])
		linha = cdsfile.readline()
	if(not found):
		logfile = open("cleaninglog.txt", "a")
		logfile.write(((splitmRNAData[2].split("="))[1].split(","))[0] + "\t" + (splitmRNAData[0].split("="))[1] + "\t" + "CDS not found\n")
		logfile.close()
	else:
		resultfile = open("mrnainfo.gff", "a")
		resultfile.write(((splitmRNAData[2].split("="))[1].split(","))[0] + "\t" + (splitmRNAData[0].split("="))[1] + "\t" + str(lowestEntry) + "\t" + str(highestEntry) + "\n")
		resultfile.close()
	linha_2 = mrnafile.readline()

print("Process Finished!")
print("Total Time Elapsed: " + str(datetime.now() - start_time))

