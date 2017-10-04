import Bio
import numpy
import argparse
from datetime import datetime
import sys
import os
import math

def MergeSort_CDS(lista):
	if(len(lista)>1):
		mid = len(lista)//2
		left = lista[:mid]
		right = lista[mid:]

		MergeSort_CDS(left)
		MergeSort_CDS(right)

		i=0
		j=0
		k=0

		while i < len(left) and j < len(right):
			if((((((left[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((left[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
				LeftID = int((((((left[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1].split(":"))[1])
			else:
				LeftID = int((((((left[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1])

			if((((((right[j].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((right[j].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
				RightID = int((((((right[j].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[1].split(":"))[1])
			else:
				RightID = int((((((right[j].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[1])

			if(LeftID < RightID):
				lista[k]=left[i]
				i=i+1
			else:
				lista[k]=right[j]
				j=j+1
			k=k+1

		while i < len(left):
			lista[k]=left[i]
			i=i+1
			k=k+1

		while j < len(right):
			lista[k]=right[j]
			j=j+1
			k=k+1
	return lista

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
parser.add_argument('-s','--sorted', dest='sort' ,action='store_true', default=False, help='mark this to skip sorting process (only if files are sorted)')


args = parser.parse_args()

start_time = datetime.now()

if(not(os.path.isfile(os.path.splitext(args.gff)[0] + "_mRNA.gff") and os.path.isfile(os.path.splitext(args.gff)[0] + "_CDS.gff"))):
	file = open(args.gff, "r")

	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA.gff", "w")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS.gff", "w")

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
else:
	print("mRNA and CDS files already exist, skipping to the next process...")

# -------------------------------------------------------

# Sorts mRNA file

if(not args.sort and not(os.path.isfile(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff"))):
	print("Starting mRNA file sorting...")
	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA.gff", "r")
	mrna = []
	linha = mrnafile.readline()

	while linha: 
		mrna.append(linha)
		linha = mrnafile.readline()

	mrnafile.close()

	mrna = QuickSort_mRNA(mrna, 0, len(mrna)-1)


	print("Writing sorted file...")
	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "w")
	for element in mrna:
		mrnafile.write(element)
	mrnafile.close()


	logfile = open("timelog.txt", "w")
	print("Time Elapsed: " + str(datetime.now() - start_time))
	logfile.write("Time Taken: " + str(datetime.now() - start_time))
else:
	print("Skipping mRNA file sorting...")

# ------------------------------------------------------
# Sorts CDS file

if(not args.sort and not(os.path.isfile(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff"))):
	print("Starting CDS file sorting...")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS.gff", "r")
	cds = []
	linha = cdsfile.readline()

	while linha: 
		cds.append(linha)
		linha = cdsfile.readline()

	cds = MergeSort_CDS(cds)

	print("Writing sorted file...")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "w")
	for element in cds:
		cdsfile.write(element)
	cdsfile.close()


	logfile = open("timelog.txt", "a")
	print("Time Elapsed: " + str(datetime.now() - start_time))
	logfile.write("Time Taken: " + str(datetime.now() - start_time + "\n"))
else:
	print("Skipping CDS file sorting...")

# ---------------------------------------------------------------

# Create partitions on the CDS file for quick indexing ----
print("Creating an entry index...")
lastelement = ''

IDlist = [0]
with open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "r") as f:
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
			#last_address = f.tell()
			#continue

		last_address = f.tell()
		linha = f.readline()
print("Size: " + str(len(IDlist)))

print("Time Elapsed: " + str(datetime.now() - start_time))

#------------------------------------------------------------
# Starts ORF extraction

mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "r")
cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "r")
resultfile = open("mrnainfo.txt", "w")
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
		resultfile = open("mrnainfo.txt", "a")
		resultfile.write(((splitmRNAData[2].split("="))[1].split(","))[0] + "\t" + (splitmRNAData[0].split("="))[1] + "\t" + str(lowestEntry) + "\t" + str(highestEntry) + "\n")
		resultfile.close()
	linha_2 = mrnafile.readline()

print("Process Finished!")
print("Total Time Elapsed: " + str(datetime.now() - start_time))

