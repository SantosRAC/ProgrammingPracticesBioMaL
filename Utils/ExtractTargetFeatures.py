import argparse
from datetime import datetime
from numpy import log10
import sys
import os

# Sorts the CDS file by GeneID
def MergeSort_CDS(lista):
	CDSData_DbxrefIndexL = -1
	CDSData_DbxrefIndexR = -1
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

			indLeft = 0
			indRight = 0

			# Finds index of Dbxref field for Left
			for element in ((left[i].split("\t"))[8].split(";")): 
				if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
					CDSData_DbxrefIndexL = indLeft
				else:
					indLeft = indLeft + 1

			# Finds index of Dbxfer field for Right
			for element in ((right[j].split("\t"))[8].split(";")): 
				if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
					CDSData_DbxrefIndexR = indRight
				else:
					indRight = indRight + 1

			if((((((left[i].split("\t"))[8].split(";"))[CDSData_DbxrefIndexL].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((left[i].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
				LeftID = int((((((left[i].split("\t"))[8].split(";"))[CDSData_DbxrefIndexL].split("="))[1].split(","))[1].split(":"))[1])
			else:
				LeftID = int((((((left[i].split("\t"))[8].split(";"))[CDSData_DbxrefIndexL].split("="))[1].split(","))[0].split(":"))[1])

			if((((((right[j].split("\t"))[8].split(";"))[CDSData_DbxrefIndexR].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((right[j].split("\t"))[8].split(";"))[2].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
				RightID = int((((((right[j].split("\t"))[8].split(";"))[CDSData_DbxrefIndexR].split("="))[1].split(","))[1].split(":"))[1])
			else:
				RightID = int((((((right[j].split("\t"))[8].split(";"))[CDSData_DbxrefIndexR].split("="))[1].split(","))[0].split(":"))[1])

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

# Sorts the mRNA file by GeneID
def MergeSort_mRNA(lista):
	RNAData_DbxrefIndexR = -1
	RNAData_DbxrefIndexL = -1
	if(len(lista)>1):
		mid = len(lista)//2
		left = lista[:mid]
		right = lista[mid:]

		MergeSort_mRNA(left)
		MergeSort_mRNA(right)

		i=0
		j=0
		k=0

		while i < len(left) and j < len(right):

			indLeft = 0
			indRight = 0

			# Finds index of Dbxref field for Left
			for element in ((left[i].split("\t"))[8].split(";")): 
				if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
					RNAData_DbxrefIndexL = indLeft
				else:
					indLeft = indLeft + 1

			# Finds index of Dbxfer field for Right
			for element in ((right[j].split("\t"))[8].split(";")): 
				if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
					RNAData_DbxrefIndexR = indRight
				else:
					indRight = indRight + 1

			LeftID = int((((((left[i].split("\t"))[8].split(";"))[RNAData_DbxrefIndexL].split("="))[1].split(","))[0].split(":"))[1])

			RightID = int((((((right[j].split("\t"))[8].split(";"))[RNAData_DbxrefIndexR].split("="))[1].split(","))[0].split(":"))[1])

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

# ------ Argsparsing ------------
sys.setrecursionlimit(1000000) # Default = 10000
parser = argparse.ArgumentParser(description='Extracts UTR length from mRNA entries', add_help=True)
parser.add_argument('-g','--gff', dest='gff', metavar='inputGFF', type=str, help='GFF containing information about features in genome', required=True)

args = parser.parse_args()
#--------------------------------

start_time = datetime.now()

# Creates separate files exclusively for mRNAs and CDS'
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


if(not(os.path.isfile(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff"))):
	print("Starting mRNA file sorting...")
	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA.gff", "r")
	mrna = []  # List containing all lines in file
	linha = mrnafile.readline()

	while linha: # Appends all lines to mrna list
		mrna.append(linha)
		linha = mrnafile.readline()		
	mrnafile.close()

	mrna= MergeSort_mRNA(mrna) # Sorts the list

	# Writes back list mrna to mRNA_Sorted file
	print("Writing sorted file...")
	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "w")
	for element in mrna:
		mrnafile.write(element)
	mrnafile.close()

	# Generates a timelog file
	logfile = open("timelog.txt", "w")
	print("Time Elapsed: " + str(datetime.now() - start_time))
	logfile.write("Time Taken: " + str(datetime.now() - start_time) + ("\n"))
else:
	print("Skipping mRNA file sorting...")

# ------------------------------------------------------
# Sorts CDS file

CDSData_Indexes = []

if(not(os.path.isfile(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff"))):
	print("Starting CDS file sorting...")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS.gff", "r")
	cds = [] # List containing all lines in file
	linha = cdsfile.readline()

	while linha: # Appends all lines to cds list
		cds.append(linha)
		linha = cdsfile.readline()
	cdsfile.close()

	cds= MergeSort_CDS(cds) # Sorts the list

	# Writes back list cds to CDS_Sorted file
	print("Writing sorted file...")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "w")
	for element in cds:
		cdsfile.write(element)
	cdsfile.close()

	# Generates a timelog file
	logfile = open("timelog.txt", "a")
	print("Time Elapsed: " + str(datetime.now() - start_time))
	logfile.write("Time Taken: " + str(datetime.now() - start_time) + "\n")
else:
	print("Skipping CDS file sorting...")

# ---------------------------------------------------------------

# Create partitions on the CDS file for quick indexing ----

if(not os.path.isfile("mrnainfo.txt")):
	print("Creating an entry index...")
	lastelement = ''  # ID Variable to check wheather there was a partition change

	IDlist = [0]  # List containing all the "addresses" of the partitions in a file
	Comparelist = []
	foundField = False

	# Creates an index list for Dbxref field
	if(os.path.isfile(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff")):
		with open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "r") as f:

			linha = f.readline()

			while linha: # Appends all lines to mrna list
				i = 0
				foundField = False
				splitString = ((linha.split("\t"))[8].split(";"))

				# Creates a list of Dbxref field indexes
				for element in splitString: 
					if(foundField):
						continue
					if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
						CDSData_Indexes.append(i)
						foundField = True
					else:
						i = i + 1
				linha = f.readline()

				if(not foundField):  # If GeneID field is not found; Kills the script
					print("Following line contains no GeneID information: " + linha)
					sys.exit()

	with open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "r") as f:
		linha = f.readline()
		j = 0

		# Adds all mRNA IDs to Compare List
		while linha:
			line_element = (((((linha.split("\t"))[8].split(";"))[CDSData_Indexes[j]].split("="))[1].split(","))[0].split(":"))[1]
			j = j + 1
			if(lastelement == ''):
				lastelement = line_element
			elif(line_element != lastelement):
				Comparelist.append(line_element)
			linha = f.readline()		

	lastelement = ''

	with open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "r") as f:
		linha = f.readline()
		j = 0

		# Compares the IDs of the actual and last element's ID. If there is a change, a new partition is added
		# Only CDS' IDs that are in Compare List are added as a partition. This operation discards all CDS data that have no related mRNA.
		while linha:
			if((((((linha.split("\t"))[8].split(";"))[CDSData_Indexes[j]].split("="))[1].split(","))[0].split(":"))[0] == "CCDS" or (((((linha.split("\t"))[8].split(";"))[CDSData_Indexes[j]].split("="))[1].split(","))[0].split(":"))[0] == "Genbank"):
				line_element = (((((linha.split("\t"))[8].split(";"))[CDSData_Indexes[j]].split("="))[1].split(","))[1].split(":"))[1]
			else:
				line_element = (((((linha.split("\t"))[8].split(";"))[CDSData_Indexes[j]].split("="))[1].split(","))[0].split(":"))[1]
			if(lastelement == ''):
				lastelement = line_element
			elif(line_element != lastelement):
				if(line_element in Comparelist):
					IDlist.append(last_address)
				lastelement = line_element
				last_address = f.tell()
				continue

			j = j + 1
			last_address = f.tell()
			linha = f.readline()

	print("Time Elapsed: " + str(datetime.now() - start_time))

else:
	print("Skipping Partitioning!")


#------------------------------------------------------------
# Starts ORF extraction


if(not os.path.isfile("mrnainfo.txt")):

	RNAData_Indexes = []

	# Creates an index list for Dbxref field on mRNA file
	if(os.path.isfile(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff")):
		with open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "r") as f:

			linha = f.readline()

			while linha: # Appends all lines to mrna list
				i = 0
				foundField = False
				splitString = ((linha.split("\t"))[8].split(";"))

				# Creates a list of Dbxref field indexes
				for element in splitString: 
					if(foundField):
						continue
					if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
						RNAData_Indexes.append(i)
						foundField = True
					else:
						i = i + 1
				linha = f.readline()

				if(not foundField):  # If GeneID field is not found; Kills the script
					print("Following line contains no GeneID information: " + linha)
					sys.exit()


	# Opens files and print their header information
	mrnafile = open(os.path.splitext(args.gff)[0] + "_mRNA_Sorted.gff", "r")
	cdsfile = open(os.path.splitext(args.gff)[0] + "_CDS_Sorted.gff", "r")
	resultfile = open("mrnainfo.txt", "w")
	resultfile.write("#GeneID" + "\t" + "rna_code" + "\t" + "cds_start_position" + "\t" + "cds_final_position" + "\t" + "strand" + "\n")
	resultfile.close()
	logfile = open("cleaninglog.txt", "w")
	logfile.write("#GeneID\trna_code\terror_message\n")
	logfile.close()

	lastGeneID = ''
	linecount = 0
	IDIndex = 0

	print("Starting mRNAs' ORF extraction...")

	line_mRNA = mrnafile.readline()

	while line_mRNA:  # Iterates through the mRNA file
		highestEntry = 0
		lowestEntry = 900000000000000000000
		sys.stdout.write('\r' + "Processing mRNA #" + str(linecount+1)) # stdout prints while erasing the last line of console
		sys.stdout.flush()
		linecount = linecount + 1
		i = 0

		# Variable that represent ID field in the gff
		RNAData_IDIndex = -1

		# Splitting of file line tabs
		splitmRNA = line_mRNA.split("\t")
		splitmRNAData = splitmRNA[8].split(";")

		for element in splitmRNAData:  # Finds index of ID field
			if(element[0] == "I" and element[1] == "D"):
				RNAData_IDIndex = i
			else:
				i = i + 1		

		# Get actual GeneID and checks if still in the same partition
		# Skips file to the next partition on the next iteration
		if(lastGeneID == ''):  
			lastGeneID = int((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(","))[0].split(":"))[1])
		elif(lastGeneID != int((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(","))[0].split(":"))[1])):
			lastGeneID = int((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(","))[0].split(":"))[1])
			IDIndex = IDIndex + 1
			cdsfile.seek(IDlist[IDIndex],0)
		else:
			cdsfile.seek(IDlist[IDIndex],0)
		
		# Logfile variables
		found = False
		nocds = True

		line_CDS = cdsfile.readline()

		while line_CDS:  # Iterates through CDS partition
			# Variables that represent Parent and Dbxref fields in the gff
			CDSData_ParentIndex = -1
			CDSData_DbxrefIndex = -1
			i = 0
			j = 0

			# Splitting of lines in CDS file
			splitCDS = line_CDS.split("\t")
			splitCDSData = splitCDS[8].split(";")

			for element in splitCDSData:  # Finds indexes of Parent and Dbxref fields
				if(element[0] == "P" and element[1] == "a" and element[2] == "r"):
					CDSData_ParentIndex = i
				else:
					i = i + 1
				if(element[0] == "D" and element[1] == "b" and element[2] == "x"):
					CDSData_DbxrefIndex = j
				else:
					j = j + 1

			if(CDSData_ParentIndex == -1 or CDSData_DbxrefIndex == -1): # If no Parent or Dbxref fields found
				line_CDS = cdsfile.readline()
				continue

			# Checks if still in the same partition
			if(((splitCDSData[CDSData_DbxrefIndex].split("="))[1].split(":"))[0] == "CCDS" or ((splitCDSData[CDSData_DbxrefIndex].split("="))[1].split(":"))[0] == "Genbank"):
				CompareID = int((((splitCDSData[CDSData_DbxrefIndex].split("="))[1].split(","))[1].split(":"))[1])
			else:
				CompareID = int((((splitCDSData[CDSData_DbxrefIndex].split("="))[1].split(","))[0].split(":"))[1])
			if(CompareID != lastGeneID):  # Stops iterating if out of partition range
				break
			# Checks if GeneID and Parent-ID information matches
			if(lastGeneID == CompareID and int(((splitmRNAData[RNAData_IDIndex].split("="))[1].split("a"))[1]) == int((((splitCDSData[CDSData_ParentIndex].split("="))[1].split("a"))[1]))):
				found = True
				if(splitCDS[3].isdigit()): # Checks if start and end positions are the lowest or highest possible
					nocds = False	
					if(int(splitCDS[3]) < lowestEntry):
						lowestEntry = int(splitCDS[3])
					if(int(splitCDS[4]) > highestEntry):
						highestEntry = int(splitCDS[4])
						
			line_CDS = cdsfile.readline()

		if(nocds): # Prints an error message in logfile
			logfile = open("cleaninglog.txt", "a")
			logfile.write((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(":"))[1].split(","))[0] + "\t" + (splitmRNAData[RNAData_IDIndex].split("="))[1] + "\t" + "No CDS positions found\n")
			logfile.close()
		elif(not found): # Prints an error message in logfile
			logfile = open("cleaninglog.txt", "a")
			logfile.write((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(":"))[1].split(","))[0] + "\t" + (splitmRNAData[RNAData_IDIndex].split("="))[1] + "\t" + "Invalid Entry\n")
			logfile.close()		
		else: # Prints all extracted results from file
			resultfile = open("mrnainfo.txt", "a")
			resultfile.write((((splitmRNAData[RNAData_Indexes[linecount-1]].split("="))[1].split(":"))[1].split(","))[0] + "\t" + (splitmRNAData[RNAData_IDIndex].split("="))[1] + "\t" + str(lowestEntry) + "\t" + str(highestEntry) + "\t")
			if(splitmRNA[6] == "+"):
				resultfile.write("plus" + "\n")
			else:
				resultfile.write("minus" + "\n")
			resultfile.close()
		line_mRNA = mrnafile.readline()


	print("Time Elapsed: " + str(datetime.now() - start_time))

else:
	print("Skipping ORF Extraction!")

# ----------------------------------------------------------------

# Gene Filtering

Genefile = open("mrnainfo.txt", "r")
print("Filtering Genes...")

Genefile.readline()  # Excludes header
line_Gene = Genefile.readline()
GeneData = line_Gene.split("\t")

# Data on the highest ranked mRNA of a gene
lastGeneID =  GeneData[0]
lastLength =  int(GeneData[3]) - int(GeneData[2])  # Subtracts final and start positions
BestGene = line_Gene

ChosenGenes = []

# Iterates through file and selects the best mRNAs
while line_Gene:
	if(lastGeneID == GeneData[0]):
		if(int(GeneData[3]) - int(GeneData[2]) > lastLength):
			BestGene = line_Gene
			lastLength = int(GeneData[3]) - int(GeneData[2])
	else:
		BestGene = BestGene[:-1]
		BestGene = BestGene + "\t" + str(lastLength) + "\t" + str(log10(int(lastLength))) + "\n"
		ChosenGenes.append(BestGene)
		lastGeneID = GeneData[0]
		lastLength = int(GeneData[3]) - int(GeneData[2])
		BestGene = line_Gene

	line_Gene = Genefile.readline()
	GeneData = line_Gene.split("\t")

Genefile.close()

# Outputs into a file
Outfile = open("FilteredGenes.txt", "w")

Outfile.write("#GeneID\trna_code\tcds_start_position\tcds_final_position\tstrand\tlength\tlog10length\n")

for line in ChosenGenes:
	Outfile.write(line)

Outfile.close()

print("Process Finished!")
print("Total Time Elapsed: " + str(datetime.now() - start_time))
