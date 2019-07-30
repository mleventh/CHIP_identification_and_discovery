#!/usr/bin/env python
#encoding: utf-8

# Not function not enabled
# fixed exac filter bug
#v2 debug -bg

#import
import re
import os
import sys
import argparse
import time

# Define command line input and help menu below
parser = argparse.ArgumentParser(
		description='Take a list of variants and filter them by list of rules',
		epilog = 'Formatting for rules file:')
parser.add_argument('rules_file', metavar='rules', type=str, 
		help = 'Path to file containing filtering rules')
parser.add_argument('infile', metavar='input', type=str, 
		help = 'Path to file containing filtering rules')
parser.add_argument('-b', '--blacklist', dest='blacklist', metavar='',
        help='Path to list of genes to exclude before applying rules')
parser.add_argument('-o', '--output', dest='outfile', 
		default='output', metavar='',
        help='Path to output file (Default: output')
parser.add_argument('-bg', '--by_gene', dest='bg', 
		default=False, metavar='', action='store_const', const=True,
        help='Create a directory of variant files sorted by gene')
args = parser.parse_args()

# Move namespace to global variables
globals().update(vars(args))


# dictionary of shorthand to function, ## to be replaced with rules details
# For annovar headers
func_lookup = {
	"position"	: "position(spl[idx('Start_position') ], spl[idx('End_position') ], ## )",
	"type" 		: "typemut(spl[idx('Variant_Classification') ], ## )",
	"aa_match" 	: "aa_match(spl[idx('Protein_Change') ], ## )",
	"aa_range"	: "aa_range(spl[idx('Protein_Change') ], ## )",
	"fs" 		: "frameshift(spl[idx('Variant_Classification') ], ## )",
	"field" 	: "fieldmatch(line, ## )",
	#"ExAC_Freq"	: "exac_freq(spl[idx('ExAC_Freq')] , ##)"
}

# Dictionaries to be loaded
rules = {  
	#gene: [rule, rule, ... etc] or as string
	}
gene_outfiles = { }
blacklist_vars = []

# Main function
def main():

	#with printing meta data can't rerun on older files
	
	# Get constructors for output files
	outdir = os.path.dirname(os.path.realpath(outfile))
	selectfile = outfile  + "_selected.txt"
	selectpath = selectfile

	dropfile = outfile + "_dropped.txt"
	droppath = dropfile

	blacklistedpath = outfile + "_blacklisted.txt"
	
	norulespath = outfile + "_norules.txt"
	
	#load rules and stats on rules file
	if blacklist != None: load_blacklist(blacklist) 
	load_rules(rules_file)


	# if by gene is true make files to write to
	if (bg == True):
		gene_dir = os.path.join(outdir, "by_gene")
		if not os.path.exists(gene_dir):
		    os.makedirs(gene_dir)

		for g in rules:
			spath = os.path.join(gene_dir, g + "_selected.txt")
			dpath = os.path.join(gene_dir, g + "_dropped.txt")
			gene_outfiles[g] = (open(spath, "w"), open(dpath, "w"))


	# Open output path for all variants
	sel = open(selectpath, "w")
	metaStamp(sel)
	drp = open(droppath, "w")
	metaStamp(drp)
	
	blk = open(blacklistedpath, "w")
	metaStamp(blk)
	nog = open(norulespath, "w")
	metaStamp(nog)


	header = []
	# open file
	try:
		with open(infile) as f:		

			# make header and function to get index my header value		
			header = f.readline().rstrip().split("\t")
			idx = lambda x: header.index(x)
			
			#print header to all output files
			for fi in gene_outfiles.values() + [(sel, drp)] + [(blk, nog)]:
				fi[0].write("\t".join(header) + "\n")
				fi[1].write("\t".join(header) + "\n")

			# go through file line by line
			for line in f:

				# Split line by tabs check if variant is in blacklist
				spl = line.rstrip().split("\t")
				var = ( spl[idx('Chromosome')], spl[idx('Start_position')], 
						spl[idx('End_position')], spl[idx('Reference_Allele')], spl[idx('Tumor_Seq_Allele2')]
					)
				if var in blacklist_vars:
					blk.write(line)
					continue

				try:
					
					# If line matches gene rules write line to selected file
					r = rules[spl[idx('Hugo_Symbol')]]
					if eval(r): 
						sel.write(line)
						if bg == True:
							gene_outfiles[spl[idx('Hugo_Symbol')]][0].write(line)

					# Else they don't match print to discarded file
					else:
						drp.write(line)
						# Print non-rules matching variables
						if bg == True:
							gene_outfiles[spl[idx('Hugo_Symbol')]][1].write(line)

				except:
					nog.write(line)
					pass

	except IOError:
		print "Could not open rule file " + rules_file

	sel.close()
	drp.close()
	blk.close()
	nog.close()

	if bg == True:
		for g in rules:
			gene_outfiles[g][0].close()
			gene_outfiles[g][1].close()



#########################################################################

# Function convert shorthand to function from lookup 
#	(inputs form csv string  eg: type:x,y,z)
def lookup (string):

	# parse by word (and parens) 
	x = re.findall(r'[\w]+[(:|,)[\w|\d|\.]+[^\)|\s]]*|\(|\)|\w+', string)

	# Ignoring parens and and/or's replace words 
	#	with correpsonging functions, populated
	for i, w in enumerate(x):
		if w != "and" and w != "or" and w != "(" and w != ")":
		# if w != "and" and w != "or" and w != "(" and w != ")" and w != "not":  # NOT STATEMENT BROKEN, needs a werid space to work
			var = w.split(':')
			x[i] = func_lookup[var[0]].replace(
					"##", "\"" + "\" , \"".join(var[1].split(",")) + "\""
				)

	# return rejoined rule
	return " ".join(x)


# Put header metadata in file
def metaStamp(f):
	rules_stats = os.stat(rules_file)
	f.write('# Filtered on %s \n' % time.strftime("%m/%d/%Y"))
	f.write('# Using rules %s \n' % rules_file)
	f.write('# Last Modified: %s \t Filesize: %s \n' % 
				(time.ctime(rules_stats.st_mtime), rules_stats.st_size)
			)

	blklist_stats = os.stat(blacklist)
	f.write('# Using blacklist %s \n' % blacklist)
	f.write('# Last Modified: %s \t Filesize: %s \n' % 
				(time.ctime(blklist_stats.st_mtime), blklist_stats.st_size)
			)


# function to load rules from file
def load_rules( rulespath ):

	# Load rules line by line
	try:
		with open(rulespath) as f:
			gene = ""
			for line in f:
				line = line.rstrip()

				# skip commented and blank lines
				if line == "" or line[0] == "#":
					continue

				# Set gene, REGEX TO JUST TAKE WORD?
				if line[0] == ">":
					gene = line[1:]
					continue

				# If gene list isn't empty add, else make and add
				try:
					rules[gene].append(lookup(line))
				except:
					rules[gene] = [lookup(line)]

	except IOError:
		print "Could not open rule file " + rulespath

	# compress rules into one line
	for g in rules:
		rules[g] = "(" + ") or (".join(rules[g]) + ")"
		print rules[g]

# function to load a blacklist from a file into the global blacklist
def load_blacklist(listpath):

	# Load list line by line
	try:
		with open(listpath) as f:
			# than put into tuples add to blacklist_vars
			for line in f:
				line = line.rstrip()
				# skip commented and blank lines
				if line == "" or line[0] == "#":
					continue

				this_var = tuple(line.split("\t"))
				blacklist_vars.append(this_var)

	except IOError:
		print "Could not open rule file " + listpath


# True/ False returning functions, on line data
#########################################################################

# given start and stop of line variant, 
#	return true if variant overlaps with given range
def position(pos1, pos2, start, finish = -1):
	pos1 = int(pos1)
	pos2 = int(pos2)
	start = int(start)
	finish = int(finish)
	if finish < 0: finish = start
	return  pos2 >= start and finish >= pos1 

# If variant  has a type matching input 
def typemut(var, typ) :
	if re.search(typ, var): return True 
	return False

# Amino acid contains given pattern treutrn true
def aa_match(line, match_string) :
	if re.search(match_string, line): return True
	else: return False

# Is there a frameshift, if yes return true
def frameshift(line, b) :
	if re.search(r"^(?!non)(?=frameshift).*", line):
		if b == "T": return True
	else:
		if b == "F": return False
	return False

# If the Amino acid in a given range/ region
def aa_range(val,start,end):
	#pos = int(re.findall(r'^p.\w(\d+)\w', val)[0])   # DOES THIS WORK FOR IN FRAME INDELS FROMAT p.360_361del try   re.findall(r'\d+', val)[0]
	pos = int(re.findall(r'\d+', val)[0]) 
	(start, end) = (int(start), int(end))
	if pos < end and pos > start:
		return True
	return False

#def exac_freq(line, thresh):
	#if float(line) < float(thresh):
		#return True
	#else:
		#return False

# def fieldmatch(line, field, val):





# Run Main function after code implementation
if __name__ == '__main__':
	main()
