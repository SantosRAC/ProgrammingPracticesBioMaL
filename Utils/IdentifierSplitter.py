# Removes all duplicated identifiers from a mirTarBase .csv and splits the output file
__author__      = "Henrique Frajacomo"

import csv      # Import .csv files
import argparse # Package for using command-line arguments in Python

parser = argparse.ArgumentParser(description='Removes all lines in which identifiers match the argument entry', add_help=True)
parser.add_argument('-f','--file', dest='file1', metavar='inputFile', type=str, help='File to be splitted', required=True)
parser.add_argument('-o','--outfilename', dest='outName', metavar='outputfileName', type=str, help='name of the output file', required=True)
parser.add_argument('-n','--number', dest='linePerFile', metavar='linePerFile', type=int, help='number of lines in each split file', required=True)

args = parser.parse_args()
file1 = args.file1
filename = args.outName
lineNumber = args.linePerFile

'''

	Functions
	
'''

# Extracts all identifiers from a .csv file and returns a list of them
def CsvExtractIdentifiers(myFile):
    Identifiers = []
    temp = ""
    for line in myFile:
        temp = line.split(";")
        Identifiers.append(temp[1])
    Identifiers.pop(0)

    return Identifiers


# Pops all duplicated elements in a list and returns the resulting list
def RemoveDuplicates(myList):
    Removed = 0
    nonDuplicatedList = []
    for element in myList:
        found = False
        ReadData = element
        for element in nonDuplicatedList:
            if(ReadData == element):
                found = True
                Removed = Removed + 1
                break
        if(not found):
            nonDuplicatedList.append(ReadData)

    print("Removed " + str(Removed) + " entries")

    return nonDuplicatedList


# Writes the content of a list to a set number of .csv files
def CsvSplitWriter(linesPerFile, baseName, myList):
    filecounter = 1
    counter = 1
    TotalLines = len(myList)

    numberOfFiles = int(TotalLines/linesPerFile)+1

    outFile = open(baseName + str(filecounter) + ".csv", "w")
    for element in myList:
        outFile.write(element + "\r")
        if(counter == linesPerFile):
            counter = 0
            filecounter = filecounter + 1
            outFile.close()
            outFile = open(baseName + str(filecounter) + ".csv", "w")
        counter = counter +1
    outFile.close()


''' 

	Main
	
'''

file1OBJ = open(file1, 'r')  # first argument .csv file

list1 = CsvExtractIdentifiers(file1OBJ)  # list containing all identifiers from the file

list1 = RemoveDuplicates(list1)  # identifiers filtered list 

CsvSplitWriter(lineNumber, filename, list1)  # Splits the list into a set number of files 

file1OBJ.close()
