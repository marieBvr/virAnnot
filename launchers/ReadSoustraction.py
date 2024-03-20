# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log


class ReadSoustraction:
	"""
	This class is a part of the virAnnot module
	It substracts reference sequences from reads
	Reference sequences must be in Bowtie format
	"""

	def __init__(self, args):
		self.execution = 1 # init value
		self.check_args(args)
		self._create_cmd()


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		log.info("Checking parameters...")
		self.wd = os.getcwd()
		self.params = args['params']
		self.cmd = []
		self.execution = 1
		if 'iter' in args:
			if args['iter'] == 'sample':
				self.sample = args['sample']
				self.cmd_file = self.wd + '/' + 'readSoustraction' + self.sample + '_cmd.txt'
				self.iter = 'sample'
			elif args['iter'] == 'library':
				self.library = args['library']
				self.cmd_file = self.wd + '/'  + 'readSoustraction' + self.library + '_cmd.txt'
				self.iter = 'library'
		if 'i1' in args:
			self.i1 = self._check_file(self.wd + '/' + args['i1'])
		else:
			log.critical('Need r1 file.')
			self.execution = 0
		if 'i2' in args:
			self.i2 = self._check_file(self.wd + '/' + args['i2'])
		else:
			log.critical('Need r2 file.')
			self.execution = 0
		if 'db' in args:
			self.db = self._check_bowtie_db(args['db'])
		else:
			log.critical('Need bowtie2 db option.')
			self.execution = 0
		if 'o1' in args:
			self.o1 = args['o1']
		else:
			log.critical('Need o1 file.')
		if 'o2' in args:
			self.o2 = args['o2']
		else:
			log.critical('Need o2 file.')
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			log.debug('n_cpu option not found. default 1')
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False


	def _create_cmd(self):
		log.info("Creating bash command...")
		cmd = '#!/bin/bash\n'
		cmd += self.params['bin']['bowtie'] + ' -p ' + str(self.n_cpu)
		cmd += ' -x ' + self.db
		cmd += ' -1 ' + self.i1 + ' -2 ' + self.i2 + ' | '
		cmd += self.params['bin']['samtools'] + ' view -bS - '
		cmd += ' | ' + self.params['bin']['samtools'] + ' view -u -f 12 -F 256 - | ' + self.params['bin']['bedtools'] + ' bamtofastq -i - -fq ' + self.o1 + ' -fq2 ' + self.o2
		log.debug(cmd)
		self.cmd.append(cmd)


	def _check_file(self, input_file):
		"""
		Verify that file exists
		"""
		try:
			open(input_file)
			return input_file
		except IOError:
			print('File not found ' + input_file)
			self.execution = 0


	def _check_bowtie_db(self, db):
		if db in self.params['ReadSoustraction']['db']:
			return self.params['ReadSoustraction']['db'][db]
		else:
			log.critical('Bowtie database not found in parameter.yaml')
