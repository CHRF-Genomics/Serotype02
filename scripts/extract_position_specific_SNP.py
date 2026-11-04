#!/bin/env python2.7
from optparse import OptionParser
from Bio import SeqIO
def main():
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option("-r", "--ref", action="store", dest="ref", help="reference (fasta)", default="")
	parser.add_option("-s", "--snps", action="store", dest="snps", help="snp table (CSV)", default="")
	parser.add_option("-v", "--vcf", action="store", dest="vcf", help="vcf file", default="")
	parser.add_option("-i", "--id", action="store", dest="sample", help="Sample ID", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file", default="Output.fasta")
	return parser.parse_args()
(options, args) = main()

output = open(options.output, 'w')
idd = str(options.sample)

# Read reference sequence
ref = SeqIO.read(options.ref, "fasta")

# Read in SNP info
f = open(options.vcf,"r")
loc = []
alt = []
for line in f:
	if "#" not in line:
		fields = line.rstrip().split("\t")
		loc.append(fields[1])
		alt.append(fields[4])

# Read in SNP locations
p = open(options.snps,"r")
snplist = []
for line in p:
	snplist.append(line.rstrip())

# Find SNP and generate sequence
seq=[]
for i in range(0,len(snplist)):
	if snplist[i] in loc:
		j = loc.index(snplist[i])
		seq.append(alt[j])
	else:
		seq.append(ref.seq[int(snplist[i])-1])

print ("Sequence Length: "+str(len(seq)))
record = '>'+idd+'\n'+(''.join(seq))
output.write(record)


