import matplotlib.pyplot as plt
import csv
import argparse
from datetime import datetime
import sys

# Radix Sort (LSD) is the fastest sorting algorithm known.
# Can only be used for integers sorting though
# Sorts an integers array. O(n) complexity.
def radixsort( aList ):
	RADIX = 10
	maxLength = False
	tmp , placement = -1, 1
	buckets = []
 
	while not maxLength:
		maxLength = True
		# declare and initialize buckets
		buckets = [list() for _ in range( RADIX )]
 
		# split aList between lists
		for	 i in aList:
			tmp = int(i / placement)
			buckets[tmp % RADIX].append( i )
			if maxLength and tmp > 0:
				maxLength = False
	 
		# empty lists into aList array
		a = 0
		for b in range( RADIX ):
			buck = buckets[b]
			for i in buck:
				aList[a] = i
				a = a + 1

	 
		# move to next digit
		placement = placement * RADIX

# ------ Argsparsing ------------
parser = argparse.ArgumentParser(description='Plots a graph of the miRNA (X-axis) and their respective amount of interations with different genes (Y-axis)', add_help=True)
parser.add_argument('-t','--tabular', dest='tab', metavar='inputTab', type=str, help='Tabular file containing all miRNAs and targets interactions', required=True)

args = parser.parse_args()
#--------------------------------

start_time = datetime.now()
file = open(args.tab, "r")
skip = True
i = 0

miRNA = []  # List of all miRNAs
Interactions = []

# Inserts all miRNA into the miRNA array
for line in file:
	i = i + 1
	sys.stdout.write('\r' + "Reading Line: " + str(i)) # stdout prints while erasing the last line of console
	sys.stdout.flush()
	if(skip):
		skip = False
		continue
	ReadData = line.split(";")[1]
	if(ReadData not in miRNA):
		miRNA.append(ReadData)

file.seek(0)
i = 0
Index = 0
skip = True
print("\n")

# Initializes interaction array with empty sets
for i in range(0,len(miRNA)):
	Interactions.append([])

# Maps all interactions
for line in file:
	i = i + 1
	Index = 0
	sys.stdout.write('\r' + "Processing Line: " + str(i)) # stdout prints while erasing the last line of console
	sys.stdout.flush()
	if(skip):
		skip = False
		continue
	ReadmiRNA = line.split(";")[1]
	ReadData = line.split(";")[4]

	Index = miRNA.index(ReadmiRNA)
	if(ReadData not in Interactions[Index]):
		Interactions[Index].append(ReadData)

i = 0
print("\n")

# Transforms all interactions into numbers
for i in range(0, len(Interactions)):
	sys.stdout.write('\r' + "Transforming Index: " + str(i)) # stdout prints while erasing the last line of console
	sys.stdout.flush()
	size = len(Interactions[i])
	if(size == 1):
		print("\n")
		print(miRNA[i])
	Interactions.pop(i)	
	Interactions.insert(i,size)
	miRNA.pop(i)	
	miRNA.insert(i,i)

print("\nTime Elapsed: " + str(datetime.now() - start_time))

radixsort(Interactions)
n = min(Interactions)
m = max(Interactions)

plt.plot(miRNA, Interactions)
plt.xlabel("miRNA")
plt.axis([0,len(miRNA), n, m])
plt.ylabel("# of interactions")
print("Maximum Interactions: " + str(max(Interactions)))
print("Minimum Interactions: " + str(min(Interactions)))
plt.show()
