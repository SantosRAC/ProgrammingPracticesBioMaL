from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Iterator = SeqIO.parse('/home/renato/tmp/tmpMeetingSaoCarlos_08042017/Data/205_1.fasta', "fasta")

table = 1
min_pro_len = 100

for record in Iterator:
	for frame in range(3):
		count=1
		for pro in record.seq[frame:].translate(table).split("*"):
			if len(pro) >= min_pro_len:
				#print ("%s...%s - length %i, frame %i" % (pro[:30], pro[-3:], len(pro), frame))
				print ("Sequence #" + str(count) + ":" + pro)
			count=count+1
