#!/usr/bin/env python
#encoding: utf-8

# Load list of preferred ref seqs, correct annotation data to correct refseq

#import
import re
import sys
import argparse

parser = argparse.ArgumentParser(
		description='Taking a variant file and a list of preferred \
			refseq gene/pairs, chancge annotation to the preferred refseq')
parser.add_argument('vars_file', metavar='variants', type=str, 
		help = 'Path to variant file')
parser.add_argument('prefRef', metavar='preferred_refseq_file', type=str, 
		help = 'Path to file containing the preferred refseqs')
parser.add_argument('-o', '--output', dest='out_file', 
		default='output.txt', metavar='',
        help='Path to output file (Default: output.txt')
args = parser.parse_args()

# Move namespace to global variables
globals().update(vars(args))


refseq_lookup = {}

# Main function
def main():

	# Load rules into a GENE -> NM_refseq dictionary
	loadPrefRefSeq(prefRef)

	# Load file to write to 
	try:
		out = open(out_file, "w")
	except IOError:
		print "Could not open %s to write to" % out_file

	# read file line by line, parse line, use gene+dict to get preferred refseq,
	#    parse annotation to get infor for preferred refseq reprint line
	try:
		with open(vars_file) as f:

			# make header and function to get index my header value		
			header = f.readline().rstrip().split("\t")
			out.write("\t".join(header) + "\n")
			idx = lambda x: header.index(x)

			# go through file line by line
			for line in f:
				spl = line.rstrip().split("\t")

				# For non-splicing lines get refseq dictionary
				if spl[idx('Variant_Classification')] != 'splicing' and spl[idx('Variant_Classification')] !=  'exonic;splicing':
					
					refseqs =  parseRefSeq(spl[idx('annotation')])

					# use gene column to look up preferred ref
					try:
						pref = refseq_lookup[spl[idx('gene')]]
						
						spl[idx('cdna')] = refseqs[pref][3]
						spl[idx('aa')] = refseqs[pref][4]
						
						new_line = "\t".join(spl) 
						out.write(new_line + "\n")
			
					except:
						out.write(line)

				# Deal with splicing mutations
				else:
					# MAKE SURE TO PRINT SPLICING BACK TO FILE
					out.write(line)

	except IOError:
		print "Could not open rule file " + sys.argv[1]

	out.close()

# Load Gene -> RefSeq from file into dictionary
def loadPrefRefSeq(path):

	# Load rules line by line
	try:
		with open(path) as f:
			
			for line in f:
				line = line.rstrip()

				# skip commented and blank lines
				if line == "" or line[0] == "#":
					continue

				try:
					[gene, pref] = line.split("\t")	
				except ValueError:
					print "the line containing '%s' is formatted wrong" % line
					continue

				refseq_lookup[gene] = pref

	except IOError:
		print "Could not open rule file " + path


# Return list of refseq tuples
def parseRefSeq(anno):
	refs = anno.split(',')
	t = {}

	refs = filter(lambda x: x != '', refs)
	
	for r in refs:
		try:
			tmp = r.split(':')
			t[tmp[1]] = tuple(tmp)
		except:
			print anno
			
	return t


# Run Main function after code implementation
if __name__ == '__main__':
	main()
