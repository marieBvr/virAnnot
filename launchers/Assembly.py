# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys
from Bio import SeqIO

class Assembly:
	"""
	Assembly step of virAnnot module
	From fasta files, assemble sequences to generate scaffolds
	Author: Sebastien Thiel
	"""

	def __init__(self, args):
		self.i1 = ""
		self.i2 = ""
		self.ising = ""
		self.execution = 1
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	# Create command
	def create_cmd(self):
		if self.prog == 'idba':
			self._create_idba_cmd()
		elif self.prog == 'spades':
			self.create_spades_cmd()
		elif self.prog == 'newbler':
			self.create_newbler_cmd()


	# Run assembly with newbler for single-end
	def create_newbler_cmd(self):
		if self.ising != '':
			cmd = '#!/bin/bash\n'
			cmd += 'runAssembly -o ' + self.wd + '/' + self.sample + '_newbler '  + self.ising 
			log.debug(cmd)
			self.cmd.append(cmd)
			cmd = 'sed -i \'s,^>contig\([0-9]*\)  length=[0-9]*   numreads=[0-9]*,>' + self.sample + '_\\1,\' ' + self.wd + '/' + self.sample + '_newbler' + '/454AllContigs.fna'
			log.debug(cmd)
			self.cmd.append(cmd)
			cmd = 'cp ' + self.wd + '/' + self.sample + '_newbler' + '/454AllContigs.fna' + ' ' + self.wd + '/' + self.out
			log.debug(cmd)
			self.cmd.append(cmd)
			self.get_singletons()
			cmd = 'cp ' + self.wd + '/' + self.sample + '_newbler/' + self.sample + '_newbler.singletons.fa' + ' ' + self.wd + '/'
			log.debug(cmd)
			self.cmd.append(cmd)
		else:
			log.critical('Must specify one ising file.')
			sys.exit(1)


	# Run assembly with Spades
	def create_spades_cmd(self):
		self.i1 = self.check_fastq_format(self.i1)
		self.i2 = self.check_fastq_format(self.i2)
		if self.ising != '':
			self.ising = self.check_fastq_format(self.ising)
		cmd = '#!/bin/bash\n'
		cmd += 'metaspades.py' + ' -1 ' + self.i1 + ' -2 ' + self.i2 + ' -o ' + self.wd + '/' + self.sample + '_spades'
		if self.ising != '':
			cmd += ' -s ' + self.ising
		log.debug(cmd)
		self.cmd.append(cmd)

		cmd = 'sed -i \'s,^>NODE_\([0-9]*\)_length_[0-9]*_cov_.*,>' + self.sample + '_\\1,\' ' + self.wd + '/' + self.sample + '_spades' + '/scaffolds.fasta'
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = 'cp ' + self.wd + '/' + self.sample + '_spades' + '/scaffolds.fasta' + ' ' + self.wd + '/' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	# Run assembly with idba
	def _create_idba_cmd(self):
		self.i1 = self.check_fastq_format(self.i1)
		self.i2 = self.check_fastq_format(self.i2)
		merged_read = os.path.basename(self.i1).split('_')[0] + '_merged.fa'
		cmd = '#!/bin/bash\n'
		cmd += 'fq2fa --merge ' + self.i1 + ' ' + self.i2 + ' ' + merged_read
		log.debug(cmd)
		self.cmd.append(cmd)
		if self.ising != '':
			self.ising = self.check_fasta_format(self.ising)
		cmd = 'idba_ud -r ' + merged_read + ' --num_threads ' + self.n_cpu + ' -o ' + self.wd + '/' + self.sample + '_idba'
		if self.ising != '':
			cmd += ' -l ' + self.ising

		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = 'sed -i \'s,^>scaffold_\([0-9]*\),>' + self.sample + '_\\1,\' ' + self.wd + '/' + self.sample + '_idba' + '/scaffold.fa'
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = 'cp ' + self.wd + '/' + self.sample + '_idba' + '/scaffold.fa' + ' ' + self.wd + '/' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	# Get singletons
	def get_singletons(self):
		# get ID of singletons
		stats_file = self.wd + '/' + self.sample + '_newbler/454ReadStatus.txt'
		f = open(stats_file, "r")
		singletons = []
		for line in f:
			kind = line.split("\t")[1].replace('\n', '')
			if kind == 'Singleton':
				singletons.append(line.split("\t")[0])
		f.close()

		# retireve sequence 
		fasta_seq = SeqIO.parse(open(self.ising), 'fasta')
		singletons_file = self.wd + '/' + self.sample + '_newbler/' + self.sample + '_newbler.singletons.fa'
		with open(singletons_file, "w") as out_file:
			for seq in fasta_seq:
				name, sequence = seq.id, str(seq.seq)
				out_file.write('>' + name + '\n')
				out_file.write(sequence + '\n')



	# Input ising format should be fasta
	# If fastq or fq, convert to fasta
	def check_fasta_format(self, in_file):
		in_file = str(object=in_file)
		self._check_file(in_file)
		out = ''
		if in_file.lower().endswith('.fq') or in_file.lower().endswith('.fastq'):
			out = os.path.splitext(in_file)[0] + '.fa'
			if not self._check_file(out):
				cmd = 'fq2fa' + ' ' + in_file + ' > ' + out
				log.debug(str(cmd))
				self.cmd.append(cmd)
			return out
		else:
			log.debug('Format seems to be fastq (ising).')
			return in_file

	# Input format (i1 and i2) should be fastq
	# If fasta, fas or fa, convert to fastq
	def check_fastq_format(self, in_file):
		in_file = str(object=in_file)
		self._check_file(in_file)
		out = ''
		if in_file.lower().endswith('.fa') or in_file.lower().endswith('.fasta') or in_file.lower().endswith('.fas'):
			out = os.path.splitext(in_file)[0] + '.fq'
			if not self._check_file(out):
				cmd = 'fasta_to_fastq' + ' ' + in_file + ' > ' + out
				log.debug(str(cmd))
				self.cmd.append(cmd)
			return out
		else:
			log.debug('Format seems to be fastq (r1 and/or r2).')
			return in_file

	# Verify that all mandatory parameters are present
	def check_args(self, args:dict):
		self.execution = 1
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_idba_cmd.txt'
		if 'i1' in args:
			self.i1 = self._check_file(self.wd + '/' + args['i1'])
		else:
			self.i1 = ''
		if 'i2' in args:
			self.i2 = self._check_file(self.wd + '/' + args['i2'])
		else:
			self.i2 = ''
		if 'ising' in args:
			self.ising = self._check_file(self.wd + '/' + args['ising'])
		else:
			self.ising = ''
		if self.i1 == '' and self.i2 == '' and self.ising == '':
			log.critical('At least one read file must be defined.')
			self.execution = 0
		if 'prog' in args:
			if args['prog'] == 'idba' or args['prog'] == 'spades' or args['prog'] == 'newbler':
				self.prog = args['prog']
			else:
				log.critical('Wrong assembly program name.')
				sys.exit(1)
		else:
			log.critical('Program name is mandatory.')
			sys.exit(1)
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = args['out']


	# Launch assembly with sge or not
	def launch(self):
		if self.sge:
			fw = open(self.cmd_file, mode='w')
			for elem in self.cmd:
				fw.write(elem + "\n")
			fw.close()
			qsub_call = "qsub -wd " + self.wd + " -V -N " + self.sample + '_idba' + ' -pe multithread ' + self.n_cpu + ' ' + self.cmd_file
			log.debug(qsub_call)
			os.system(qsub_call)
		else:
			for elem in self.cmd:
				os.system(elem)


	# Exisiting file
	def _check_file(self, input_file):
		try:
			open(input_file)
			return input_file
		except IOError:
			print('File not found ' + input_file)
			self.execution = 0
