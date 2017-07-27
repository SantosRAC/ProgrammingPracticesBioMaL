import Bio
import argparse
import re



parser = argparse.ArgumentParser(description='Extract target annotation from alignment')
parser.add_argument('-b','--blastResult', dest='blastResult', metavar='Result.blast', type=str, help='Tabular file (6 std qlen slen) with BLAST (BLAST 2.6.0+) results', required=True)
parser.add_argument('-gff', dest='gff', metavar='file.gff', type=str, help='Annotation file (GFF3)', required=True)
args = parser.parse_args() # Creating a list ('args') with processed arguments

blastResultOBJ = args.blastResult # Creating a python object representing BLAST results, imported using the package argparse
gffOBJ = args.gff

blastres=open(blastResultOBJ)

for line in blastres:
 [qaccver, saccver, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, sstrand]=line.split('\t')
 

blastres.close()

gfffile=open(gffOBJ)

for line in gfffile:
 p = re.compile('#')
 if p.match(line):
  #print(line)
  continue
 else:
  [idseq, source, ftType, pStart, pEnd, score, strand, ph, att]=line.split('\t')
  p=re.compile('GeneID:(\d+)')
  if ftType == 'gene':
   idGene=p.findall(att).pop()
   print("Gene: ",idGene)
  elif ftType == 'CDS':
   idGene=p.findall(att).pop()
   print("CDS in gene: ",idGene)
  elif ftType == 'transcript':
   idGene=p.findall(att).pop()
   print("Transcript in gene: ",idGene)
  else:
   continue
  
gfffile.close()
