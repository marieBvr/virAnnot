"""
Authors: Sebastien Theil, Marie Lefebvre
"""
# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import sys


class Rps2tree:
	"""
	This class is a part of the virAnnot module.
	It creates the command that will run rps perl module.
	The perl module generate OTUs from Blastx ad Blastn results.
	"""

	def __init__(self, args):
		self.execution = 1
		self.check_args(args)
		self.cmd = []
		self.create_cmd()


	def create_cmd(self):
		"""
		Create command
		"""
		cmd = '#!/bin/sh\n'
		cmd += 'rps2tree.py'
		if self.iter == 'global':
			for s_id in self.blast_files:
				cmd += ' -f ' + self.blast_files[s_id]['contigs']
				cmd += ' -r ' + self.blast_files[s_id]['pfam']
				cmd += ' -b ' + self.blast_files[s_id]['blast']
				cmd += ' -c ' + self.blast_files[s_id]['rn']
		else:
			log.debug('msg')
			sys.exit(0)
		cmd += ' -mpl ' + str(self.min_prot)
		cmd += ' -vp ' + str(self.viral_portion)
		cmd += ' -p ' + str(self.perc)
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		self.wd = os.getcwd()
		self.params = args['params']
		self.cmd_file = self.wd + '/' + 'rps2tree_cmd.txt'
		if 'out' in args:
			self.out = args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'viral_portion' in args:
			self.viral_portion = args['viral_portion']
		if 'perc' in args:
			self.perc = args['perc']
		if 'min_prot' in args:
			self.min_prot = args['min_prot']
		if 'iter' in args:
			if args['iter'] == 'global':
				self.iter = 'global'
				self.blast_files = {}
				for s_id in args['args']:
					if s_id not in self.blast_files:
						if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']) and os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['blast']) and os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id]['contigs']):
							self.blast_files[s_id] = {}
							self.blast_files[s_id]['pfam'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['pfam']
							self.blast_files[s_id]['blast'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['blast']
							self.blast_files[s_id]['contigs'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['contigs']
							self.blast_files[s_id]['rn'] = self.wd + '/' + s_id + '/' + args['args'][s_id]['rn']
		else:
			log.critical('No iter parameters.')
		if len(self.blast_files.keys()) == 0:
			self.execution = 0
