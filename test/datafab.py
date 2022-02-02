import argparse
import numpy
import random

parser = argparse.ArgumentParser(
	description='Data fabricator for testing rocketchip')
parser.add_argument('name', type=str, metavar='<name>',
	help='base name of files to create (*.fa *.fastq)')
parser.add_argument('stacks', type=int, metavar='<stacks>',
	help='number of read stacks')
parser.add_argument('rps', type=int, metavar='<rps>',
	help='number of reads per stack')
parser.add_argument('--padding', type=int, metavar='<int>', required=False,
	default=1000, help='spacing between read stacks [%(default)i]')
parser.add_argument('--stddev', type=float, metavar='<float>', required=False,
	default=0.1, help='scatter of reads [%(default).3f]')
parser.add_argument('--width', type=int, metavar='<int>', required=False,
	default=0, help='width of read stack [%(default)i]')
parser.add_argument('--length', type=int, metavar='<int>', required=False,
	default=100, help='length of reads [%(default)i]')
parser.add_argument('--paired', type=int, metavar='<int>', required=False,
	default=0, help='generate paired reads spanning <int> nt')
parser.add_argument('--seed', type=int, metavar='<int>', required=False,
	help='use random seed')
parser.add_argument('--flank', type=int, metavar='<int>', required=False,
	default=1000, help='flanking sequence [%(default)i]')
arg = parser.parse_args()

if arg.seed:
	random.seed(arg.seed)
	numpy.random.seed(arg.seed)

###############################
## Create positions of reads ##
###############################

coor = arg.flank
locs = [] # location of each read (on center I suppose)
for i in range(arg.stacks):
	for r in numpy.random.normal(0, arg.stddev, arg.rps):
		pos = int(r*arg.length)
		if arg.width > 0: pos = random.randint(pos-arg.width, pos+arg.width)
		locs.append(coor + pos)
	coor += arg.padding
coor += arg.flank

###################
## Create genome ##
###################

genome = ''
with open(f'{arg.name}.fa', 'w') as fp:
	fp.write(f'>{arg.name}\n')
	count = 0
	line = 80
	for i in range(coor):
		count += 1
		nt = random.choice('ACGT')
		fp.write(nt)
		if count % line == 0: fp.write('\n')
		genome += nt

##################
## Create fastq ##
##################

def write_fastq(fp, n, seq):
	fp.write(f'@read.{n}\n{seq}\n+\n')
	fp.write('F' * len(seq))
	fp.write('\n')

half = arg.length // 2 # half single
n = 0

if arg.paired:
	fp1 = open(f'{arg.name}_1.fastq', 'w')
	fp2 = open(f'{arg.name}_2.fastq', 'w')
else:
	fp = open(f'{arg.name}.fastq', 'w')

for pos in locs:
	n += 1
	if arg.paired:
		p1 = pos - arg.paired // 2
		p2 = pos + arg.paired // 2
		seq1 = genome[p1-half:p1+half]
		seq2 = genome[p2-half:p2+half]
		write_fastq(fp1, n, seq1)
		write_fastq(fp2, n, seq2)
	else:
		seq = genome[pos-half:pos+half]
		write_fastq(fp, n, seq)

