import sys
from Bio import SeqIO
from argparse import (ArgumentParser, FileType)

def parse_args():
	"Parse the input arguments, use '-h' for help"
	commands = ArgumentParser(description='Detect phage, repeat and variation regions in genome.')
	commands.add_argument('--genbank', type=str, required=True,
						help='Genome in Genbank format (Required).')
	return commands.parse_args()
args = parse_args()

key1 = 'bacteriophage'
key2 = 'phage'

for rec in SeqIO.parse(open(args.genbank,"r"), 'genbank'):
	feat = len(rec.features)

	for x in range(0, feat):

		if rec.features[x].type == 'repeat_region' or rec.features[x].type == 'variation':
			print (str(rec.features[x].location).replace(':',',').replace('[','').replace('](+)','').replace('](-)','') +'\t'+str(rec.features[x].type))

		if (rec.features[x].type == 'CDS'):
			try:
				if ((key1 in str(rec.features[x].qualifiers['product'])) or (key1 in str(rec.features[x].qualifiers['note'])) or (key2 in str(rec.features[x].qualifiers['product'])) or (key2 in str(rec.features[x].qualifiers['note']))):
					print (str(rec.features[x].location).replace(':',',').replace('[','').replace('](+)','').replace('](-)','')+'\t'+str(rec.features[x].type)+'\t'+str(rec.features[x].qualifiers['product'])+'\t'+str(rec.features[x].qualifiers['note']))

			except KeyError:
				if ((key1 in str(rec.features[x].qualifiers['product'])) or (key2 in str(rec.features[x].qualifiers['product']))):
					print (str(rec.features[x].location).replace(':',',').replace('[','').replace('](+)','').replace('](-)','')+'\t'+str(rec.features[x].type)+'\t'+str(rec.features[x].qualifiers['product']))

