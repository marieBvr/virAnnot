"""
Authors         :Sebastien Theil, Marie Lefebvre
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log


class Ecsv2krona:
	"""
	This class is part of virAnnot module
	It creates the command that will genrate
	a krona chart
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
		cmd = '#!/bin/bash\n'
		cmd += 'ecsv2krona.py'
		for s_id in self.blast_files:
			for i in range(0, len(self.blast_files[s_id]['id'])):
				cmd += ' -i ' + self.blast_files[s_id]['id'][i]
				cmd += ' ' + self.blast_files[s_id]['csv_file'][i]
		if self.out != '':
			cmd += ' -o ' + self.out
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.execution = 1
		if 'out' in args:
			self.out = args['out']
		else:
			self.out = './krona'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			self.n_cpu = '1'
		if 'iter' in args:
			if args['iter'] == 'global':
				self.wd = os.getcwd()
				self.cmd_file = self.wd + '/' + 'ecsv2krona_cmd.txt'
				if 'out' in args:
					self.out = args['out']
				self.iter = 'global'
				self.blast_files = {}
				for s_id in args['args']:
					if s_id not in self.blast_files:
						self.blast_files[s_id] = {}
						self.blast_files[s_id]['csv_file'] = []
						self.blast_files[s_id]['id'] = []
					for i in range(1, 100, 1):
						id_name = 'id' + str(object=i)
						opt_name = 'b' + str(object=i)
						if id_name not in args['args'][s_id] and opt_name not in args['args'][s_id]:
							continue
						if opt_name in args['args'][s_id]:
							if os.path.exists(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name]):
								self.blast_files[s_id]['csv_file'].append(self.wd + '/' + s_id + '/' + args['args'][s_id][opt_name])
								self.blast_files[s_id]['id'].append(args['args'][s_id][id_name])
			elif args['iter'] == 'sample':
				self.sample = args['sample']
				self.wd = os.getcwd() + '/' + self.sample
				self.cmd_file = self.wd + '/' + 'ecsv2krona_cmd.txt'
				if 'out' in args:
					self.out = args['out']
				self.iter = 'sample'
				self.blast_files = {}
				if self.sample not in self.blast_files:
					self.blast_files[self.sample] = {}
					self.blast_files[self.sample]['csv_file'] = []
					self.blast_files[self.sample]['id'] = []
				for i in range(1, 100, 1):
					id_name = 'id' + str(object=i)
					opt_name = 'b' + str(object=i)
					if id_name not in args and opt_name not in args:
						continue
					if os.path.exists(self.wd + '/' + args[opt_name]):
						if opt_name in args:
							self.blast_files[self.sample]['csv_file'].append(self._check_file(self.wd + '/' + args[opt_name]))
							self.blast_files[self.sample]['id'].append(args[id_name])
				if len(self.blast_files[self.sample]['csv_file']) == 0:
					self.execution = 0
			else:
				log.error('iter argument not handled. ' + args['iter'])
				self.execution = 0
		if len(self.blast_files.keys()) == 0:
			self.execution = 0


	def _check_file(input_file):
		"""
		Verify that file exists
		"""
		try:
			open(input_file)
			return input_file
		except IOError:
			print('File not found ' + input_file)
			self.execution = 0
