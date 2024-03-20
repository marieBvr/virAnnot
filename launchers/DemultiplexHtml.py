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

class DemultiplexHtml:
	"""
	This module is part of virAnnot module
	It creates the command that will generate
	the HTML report of Demultiplex step
	"""

	def __init__(self, args):
		self.check_args(args)
		self.cmd = []
		self._create_cmd()


	def _create_cmd(self):
		"""
		Create command
		"""
		cmd = '#!/bin/sh\n'
		cmd += 'demultiplex_html.py'
		keys = sorted(self.lib)
		for lib_name in keys:
			cmd += ' -i ' + lib_name + ' -c ' + self.lib[lib_name]
		cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args):
		"""
		Check if arguments are valid
		"""
		self.wd = os.getcwd()
		self.cmd_file = self.wd + '/' + 'demultiplexHtml_cmd.txt'
		self.execution = 1
		if 'out' in args:
			self.out = args['out']
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = args['n_cpu']
		else:
			self.n_cpu = str(1)
		if 'iter' in args:
			if args['iter'] == 'global':
				self.iter = 'global'
				self.lib = {}
				for s_id in args['args']:
					self.lib[args['args'][s_id]['id']] = args['args'][s_id]['csv']
			else:
				log.critical('iter parameter must be global.')
				self.execution = 0
		else:
			log.critical('No iter parameters.')
			self.execution = 0
