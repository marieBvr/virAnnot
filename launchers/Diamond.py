"""
This class is a part of the virAnnot module
Authors: Sebastien Theil, Marie Lefebvre
"""
# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys

class Diamond:
	"""
	This Class is part of virAnnot module
	It creates the command that will run DIAMOND
	http://dx.doi.org/10.1038/nmeth.3176
	"""

	def __init__(self, args):
		log.info('Diamond module')
		self.cmd = []
		self.check_args(args)
		self._create_cmd()

	def _create_cmd(self):
		cmd = '#!/bin/bash\n'
		merged_read = self.wd + '/' + os.path.basename(self.i1).split('_')[0] + '_mergedPair.fa'
		fasta_singleton = self.wd + '/' + os.path.basename(self.ising).split('.')[0] + '.fa'
		if not os.path.exists(merged_read):
			if self.i1 != '' and self.i2 != '':
				cmd += 'fq2fa --merge ' + self.i1 + ' ' + self.i2 + ' ' + merged_read
				log.debug(cmd)
				self.cmd.append(cmd)
			if self.ising != '':
				cmd += 'fq2fa ' + self.ising + ' ' + fasta_singleton
				log.debug(cmd)
				self.cmd.append(cmd)
			if self.i1 != '' and self.i2 != '' and self.ising != '':
				cmd += 'cat ' + fasta_singleton + ' >> ' + merged_read
				log.debug(cmd)
				self.cmd.append(cmd)
		# outfmt 5 = XML format
		cmd += 'diamond blastx --outfmt 5'
		cmd += ' --db ' + str(self.db)
		cmd += ' --out ' + self.out
		cmd += ' --threads ' + self.n_cpu
		# Path to the query input file
		if self.i1 != '' and self.i2 != '':
			cmd += ' --query ' + merged_read
		elif self.ising != '':
			cmd += ' --query ' + fasta_singleton
		elif self.contigs != '':
			cmd += ' --query ' + self.contigs
		cmd += ' --evalue ' + self.evalue
		cmd += ' --max-target-seqs ' + self.max_target_seqs
		cmd += ' --min-score ' + self.score
		# Report only alignments above the given percentage of sequence identity
		cmd += ' --id ' + self.identity
		# Report only alignments above the given percentage of query cover
		cmd += ' --query-cover ' + self.qov
		# Report only alignments above the given percentage of subject cover
		cmd += ' --subject-cover ' + self.hov
		if self.sensitive:
			cmd += ' --sensitive '
		if self.more_sensitive:
			cmd += ' --more-sensitive '
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		if 'iter' in args:
			if args['iter'] == 'library':
				self.library = args['library']
				self.wd = os.getcwd()
				self.iter = 'library'
				self.cmd_file = self.library + '_dmd_cmd.txt'
			elif args['iter'] == 'sample':
				self.sample = str(args['sample'])
				self.wd = os.getcwd() + '/' + self.sample
				self.iter = 'sample'
				self.cmd_file = self.sample + '_dmd_cmd.txt'
		else:
			self.sample = str(args['sample'])
			self.wd = os.getcwd() + '/' + self.sample
			self.iter = 'sample'
			self.cmd_file = self.sample + '_dmd_cmd.txt'
		if 'contigs' in args:
			self.contigs = args['contigs']
		else:
			self.contigs = ''
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
		if self.i1 == '' and self.i2 == '' and self.ising == '' and self.contigs == '':
			self.execution = 0
			log.critical('At least one read file must be defined.')
		if 'params' in args:
			self.params = args['params']
		else:
			sys.exit('Parameters not found.')
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			log.debug('n_cpu option not found. default 1')
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'sensitive' in args:
			self.sensitive = bool(args['sensitive'])
		else:
			self.sensitive = False
		if 'more_sensitive' in args:
			self.more_sensitive = bool(args['more_sensitive'])
		else:
			self.more_sensitive = False
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
		else:
			self.out = self.wd + '/' + self.sample + '_' + 'dmd.xml'
		if 'score' in args:
			self.score = str(args['score'])
		else:
			self.score = '0'
		if 'max_target_seqs' in args:
			self.max_target_seqs = str(args['max_target_seqs'])
		else:
			self.max_target_seqs = '5'
		if 'evalue' in args:
			self.evalue = str(args['evalue'])
		else:
			self.evalue = '0.1'
		if 'identity' in args:
			self.identity = str(args['identity'])
		else:
			self.identity = '0'
		if 'qov' in args:
			self.qov = str(args['qov'])
		else:
			self.qov = '0'
		if 'hov' in args:
			self.hov = str(args['hov'])
		else:
			self.hov = '0'
		if 'db' in args:
			if args['db'] in self.params['Diamond']['db']:
				self.db = self.params['Diamond']['db'][args['db']]
			else:
				self.execution = 0
				sys.exit('Database ' + args['db'] + ' not in parameters.yaml.')
		else:
			sys.exit('You must provide a database.')
		self.execution = 1


	def _check_file(self, file_path):
		"""
		Check that file exists
		"""
		try:
			open(file_path)
			return file_path
		except IOError:
			print('File not found ' + file_path)
			self.execution = 0
