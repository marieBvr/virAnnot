"""
	This class is a part of the virAnnot module
	===========
	Authors: Sebastien Theil, Marie Lefebvre
	Objective: Run the automapper toolkit
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys

class Automapper:

	def __init__(self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()

	def create_cmd(self):
		"""
		Generate command that will be ran
		"""
		cmd = '#!/bin/bash\n'
		cmd += 'autoMapper.pl'
		cmd += ' -r ' + self.ref
		cmd += ' -q ' + self.contigs
		cmd += ' -b ' + self.ecsv
		cmd += ' -o ' + self.out
		cmd += ' -p '
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Verify that all mandatory arguments are present
		"""
		self.execution = 1
		if 'sample' in args:
			self.sample = args['sample']
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + 'autoMapper_cmd.txt'
		if 'params' in args:
			self.params = args['params']
		else:
			sys.exit('Parameters not found.')
		if 'out' in args:
			self.out = self.wd + '/' + args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'ref' in args:
			if args['ref'] in self.params['servers']['enki']['db']:
				self.ref = self.params['servers']['enki']['db'][args['ref']]
			else:
				sys.exit('Database ' + args['ref'] + ' not in parameters.yaml.')
		else:
			sys.exit('You must provide a reference Blast database.')
		if 'contigs' in args:
			if os.path.exists(self.wd + '/' + args['contigs']):
				self.contigs = self.wd + '/' + args['contigs']
			else:
				self.contigs = ''
				self.execution = 0
		if 'ecsv' in args:
			if os.path.exists(self.wd + '/' + args['ecsv']):
				self.ecsv = self.wd + '/' + args['ecsv']
			else:
				self.ecsv = ''
				self.execution = 0


	def _check_file(filepath):
		""" 
		Check if file exists
		"""
		try:
			open(filepath)
			return filepath
		except IOError:
			print('File not found ' + filepath)


	def check_seq_format(self, in_file):
		""" 
		Input format must be fastq
		If fa, fasta or fas format, convert to fastq
		"""
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
			log.debug('Format seems to be fastq.')
			return in_file
