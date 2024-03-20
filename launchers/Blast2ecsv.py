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
		cmd += 'blast2ecsv.py'
		cmd += ' -a ' + self.type # algo
		cmd += ' -x ' + self.x # xml
		cmd += ' -c ' + self.contigs # contigs
		cmd += ' -me ' + self.evalue
		cmd += ' -o ' + self.out
		cmd += ' -rn ' + self.rn
		if self.score != '':
			cmd += ' -s ' + self.score
		if self.qov != 0:
			cmd += ' -qov ' + self.qov
		if self.hov != 0:
			cmd += ' -mhov ' + self.hov
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
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'out' in args:
			self.out = args['out']
		if 'params' in args:
			self.params = args['params']
		if 'evalue' in args:
			self.evalue = str(args['evalue'])
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'contigs' in args:
			if os.path.exists(self.wd + '/' + args['contigs']):
				self.contigs = self.wd + '/' + args['contigs']
			else:
				log.critical('Blast file does not exists')
				sys.exit(1)
		else:
			log.critical('You must provide a fasta blast file.')
			sys.exit(1)
		if 'rn' in args:
			if os.path.exists(self.wd + '/' + args['rn']):
				self.rn = self.wd + '/' + args['rn']
			else:
				log.critical('RN file does not exists')
				sys.exit(1)
		else:
			log.critical('You must provide a Read count file.')
			sys.exit(1)
		if 'xml' in args:
			if os.path.exists(self.wd + '/' + args['xml']):
				self.x = self.wd + '/' + args['xml']
			else:
				log.critical('XML file does not exists')
				sys.exit(1)
		else:
			log.critical('You must provide a XML result blast file.')
			sys.exit(1)
		if 'type' in args:
			self.type = args['type']
		else:
			log.error('You must provide a Blast type.')
			sys.exit(1)
		if 'score' in args:
			self.score = str(args['score'])
		else:
			self.score = '0'
		if 'qov' in args:
			self.qov = str(args['qov'])
		else:
			self.qov = '0'
		if 'hov' in args:
			self.hov = str(args['hov'])
		else:
			self.hov = '0'

