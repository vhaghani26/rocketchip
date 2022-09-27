import sys
import os

path = sys.argv[1]

WD = os.getcwd()
with open(path) as fh:
	while True:
		line = fh.readline()
		if line == '': break
		line = line.strip()
		if line.startswith('&'): continue
		if line.startswith('#'): 
			layout = line[1:]
			os.makedirs(layout, exist_ok=True)
		if line.startswith('test'):
			test = line.replace(' ', '_')
			DATA = f'{WD}/{layout}/{test}'
			os.makedirs(DATA, exist_ok=True)
		f = line.split()
		if len(f) == 9:
			name, stacks, rps, padding, stddev, width, length, paired, flank = map(str, f)
			os.chdir(f'{WD}/..')
			if os.path.isfile(f'{DATA}/genome.fa'):
				cmd = f'./datasynth {name} {stacks} {rps} -o {DATA} --genome {DATA}/genome.fa --padding {padding} --stddev {stddev} --width {width} --length {length} --flank {flank}'
				if layout.startswith('paired'): cmd += f' --paired {paired}'
				elif layout.startswith('single'): pass
				else:
					sys.exit(f'Layout error: {test}')
				os.system(cmd)
				sys.stderr.write(cmd + '\n')
			else:
				genome = f'{name}.fa'
				cmd = f'./datasynth {name} {stacks} {rps} -o {DATA} --control input --padding {padding} --stddev {stddev} --width {width} --length {length} --flank {flank}'
				if layout.startswith('paired'): cmd += f' --paired {paired}'
				elif layout.startswith('single'): pass
				else: sys.exit(f'Layout error: {test}')
				os.system(cmd)
				sys.stderr.write(cmd + '\n')
				os.chdir(DATA)
				os.system(f'mv {genome} genome.fa')
				sys.stderr.write(f'mv {genome} genome.fa\n')

