import Bio
import numpy
import argparse
import re
import orf_finder

parser = argparse.ArgumentParser(description='Extracts features of target sequences (GFF | FASTA)', add_help=True)
parser.add_argument('-g','--gff', dest='gff', metavar='inputGFF', type=str, help='GFF containing information about features in genome', required=True)

args = parser.parse_args()



InputGFF = open(args.gff, "r")

for line in InputGFF:
	p = re.compile('#')
 	if p.match(line):
		continue

	splitList = line.split('\t')

	if (splitList[2] == 'gene'):
		print(line)
