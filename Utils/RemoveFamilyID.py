# Removes all lines in which identifiers match the second argument entry
__author__      = "Henrique Frajacomo"

import csv # Import .csv files
import argparse # Package for using command-line arguments in Python

parser = argparse.ArgumentParser(description='Removes all lines in which identifiers match the argument entry', add_help=True)
parser.add_argument('-f','--file', dest='file1', metavar='inputFile', type=str, help='File to be cured', required=True)
parser.add_argument('-s','--string', dest='family', metavar='familyName', type=str, help='family prefix to be removed', required=True)
parser.add_argument('-o','--outfilename', dest='outName', metavar='outputfileName', type=str, help='name of the output file', required=True)

args = parser.parse_args()
file1 = args.file1
family = args.family
filename = args.outName

'''

	Functions
	
'''

# Extracts all information from a .csv file and returns a list
def CsvFullExtract(myFile):
    Identifiers = []
    for line in myFile:
        Identifiers.append(line.rstrip('\n'))

    return Identifiers
	
# Pops all elements whose identifiers match the first argument entry
def RemoveFamilyID(familyName, myList):
    i = 0
    Removed = 0
    temp = ""
    outList = []
    for element in myList:
        temp = element.split(";")
        if(familyName in temp[1]):
            Removed = Removed + 1
        else:
            outList.append(element)
            i = i + 1

    print("Removed " + str(Removed) + " entries")

    return outList


# Writes a list of identifiers to a Csv file as a list of identifiers
def CsvWrite(myList, filename):
    outFile = open(filename + ".csv", "w")
    for element in myList:
        outFile.write(element + "\r")
    outFile.close()

'''

    Main

'''


file1OBJ = open(file1, "r")  # first argument .csv file

list1 = CsvFullExtract(file1OBJ)     # list containing all the data from the file
list1 = RemoveFamilyID(family, list1)  # removes all matches from the list

CsvWrite(list1, filename)   # writes the list as a .csv file

file1OBJ.close()

