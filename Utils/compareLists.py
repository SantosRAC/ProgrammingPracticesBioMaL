# Returns a csv file of commom elements between two files
__author__      = "Henrique Frajacomo & Renato Augusto Correa dos Santos"

import csv  # Import CSV files
import argparse # Package for using command-line arguments in Python

parser = argparse.ArgumentParser(description='Compare lists of identifiers', add_help=True)
parser.add_argument('-l1','--list1', dest='list1', metavar='inputList1', type=str, help='List of identifiers 1', required=True)
parser.add_argument('-l2','--list2', dest='list2', metavar='inputList2', type=str, help='List of identifiers 2', required=True)
parser.add_argument('-o','--output', dest='list3', metavar='outputList', type=str, help='Resulting list (common elements)', required=True)

args = parser.parse_args()
list1 = args.list1
list2 = args.list2
list3 = args.list3

'''

                            Functions

'''
# Read a FASTA with microRNAs and returns a list with identifiers
def Extract_fasta_miRNA_ID(myFile, myList):
    i = True
    for line in myFile:
        for word in line.split():
            if(i):
                myList.append(word[1:].lower())
            break
        i = not i
    return myList

# Read a txt with microRNAs and returns a list with identifiers
def Extract_txt_miRNA_ID(myFile, myList):
    for line in myFile:
        myList.append(line[:-1])
    return myList

# Read microRNA identifiers in CSV file and returns a list with these IDs
def Extract_csv_miRNA_ID(myFile, myList):
    x = 0
    for line in myFile:
        if(x>0):
        	myList.append(line.strip("\n").lower())
        else:
        	x = x+1
    return myList

# Compares the lists and returns a list with common elements
def ListInList(myList1, myList2):
    listResults = []

    x = 0
    for x in range(0, len(myList2)):
        if(myList2[x].split(";")[1] in myList1):
            listResults.append(myList2[x])
        x = x + 1
    return listResults

# Takes all the elements in a list 2 out of list 1
def ListSubstraction(list1, list2):
    return set(list1) - set(list2)

# TODO
def PrintListtoFile(myList, myFile):
    x = 0
    for x in range(0, len(myList)):
        myFile.write(myList[x] + '\n')

'''

                        Main

'''

line1 = []  # All lines in input FASTA
line2 = []  # All identifiers in list (CSV) of identifiers
ResultLine = []  # Stores all common miRNA identifiers

list1OBJ=open(list1,"r")
list2OBJ=open(list2,"r")
list3OBJ=open(list3,"w")

line1 = Extract_fasta_miRNA_ID(list1OBJ, line1)
line2 = Extract_csv_miRNA_ID(list2OBJ, line2)
resultLine = ListInList(line1, line2)
#resultLine = list(set(line1)&set(line2))
PrintListtoFile(resultLine, list3OBJ)

list1OBJ.close()
list2OBJ.close()
list3OBJ.close()
