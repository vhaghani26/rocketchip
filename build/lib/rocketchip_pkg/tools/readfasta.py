import gzip
import random
import math
import sys


def read_record(filename):
	"""
	Input: path to fasta file
	Return: generates a tuple of (name, seq)
	Example:
		for name, seq in readfasta.read_record(filepath):
			print(name, seq)
	"""

	label = None
	seq = []

	fp = None
	if    filename == '-':         fp = sys.stdin
	elif filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seq) > 0:
				seq = ''.join(seq)
				yield(label, seq)
				label = line[1:]
				seq = []
			else:
				label = line[1:]
		else:
			seq.append(line)
	yield(label, ''.join(seq))
	fp.close()
