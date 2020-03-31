"""
	This class is a part of the virAnnot module
	=============
	Authors: Sebastien Theil, Marie Lefebvre
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys

class Blast2ecsv:
	"""
	This class is part of virAnnot module
	It creates the command that will convert
	Blast results xml format to csv format
	"""

	def __init__(self, args):
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command
		"""
		cmd = '#!/bin/bash\n'
		cmd += 'blast2ecsv.pl'
		cmd += ' -t ' + self.type
		cmd += ' -b ' + self.b
		cmd += ' -if ' + self.in_format
		cmd += ' -e ' + str(self.evalue)
		cmd += ' -pm ' + self.pm
		if self.fhit:
			cmd += ' -fhit '
		if self.contigs != '':
			cmd += ' -seq ' + self.contigs
		if self.vs:
			cmd += ' -vs '
		if self.r:
			cmd += ' -r '
		cmd += ' -o ' + self.out
		if self.rn != '':
			cmd += ' -rn ' + self.rn
		if self.score != '':
			cmd += ' -s ' + self.score
		if self.qov != 0:
			cmd += ' -qov ' + self.qov
		if self.hov != 0:
			cmd += ' -hov ' + self.hov
		if self.identity != 0:
			cmd += ' -identity ' + self.identity
		if self.pd:
			cmd += ' -pd '
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		if 'iter' in args:
			if args['iter'] == 'library':
				self.library = args['library']
				self.wd = os.getcwd()
				self.iter = 'library'
				self.cmd_file = self.library + '_b2e_cmd.txt'
			elif args['iter'] == 'sample':
				self.sample = str(args['sample'])
				self.wd = os.getcwd() + '/' + self.sample
				self.iter = 'sample'
				self.cmd_file = self.sample + '_b2e_cmd.txt'
		else:
			self.sample = str(args['sample'])
			self.wd = os.getcwd() + '/' + self.sample
			self.iter = 'sample'
			self.cmd_file = self.sample + '_b2e_cmd.txt'
		if 'contigs' in args:
			self.contigs = self.wd + '/' + args['contigs']
		else:
			self.contigs = ''
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = args['out']
		if 'params' in args:
			self.params = args['params']
		if 'evalue' in args:
			self.evalue = args['evalue']
		else:
			self.evalue = 0.0001
		if 'fhit' in args:
			self.fhit = bool(args['fhit'])
		else:
			self.fhit = False
		if 'fhsp' in args:
			self.fhit = bool(args['fhsp'])
		else:
			self.fhsp = False
		if 'pm' in args:
			self.pm = args['pm']
		else:
			self.pm = 'local'
		if 'if' in args:
			self.in_format = args['if']
		else:
			self.in_format = 'xml'
		if 'r' in args:
			self.r = bool(args['r'])
		else:
			self.r = False
		if 'vs' in args:
			self.vs = bool(args['vs'])
		else:
			self.vs = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'b' in args:
			if os.path.exists(self.wd + '/' + args['b']):
				self.b = self.wd + '/' + args['b']
			else:
				self.b = ''
				self.execution = 0
		else:
			log.critical('You must provide a blast file.')
			sys.exit(1)
		if 'rn' in args:
			self.rn = self.wd + '/' + args['rn']
		else:
			self.rn = ''
		if 'type' in args:
			self.type = args['type']
		else:
			log.error('You must provide a Blast type.')
			sys.exit(1)
		if 'score' in args:
			self.score = str(args['score'])
		else:
			self.score = '0'
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
		self.pd = 'pd' in args
