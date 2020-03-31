"""
	This class is a part of the virAnnot module
	=========
	Authors: Sebastien Theil, Marie Lefebvre
"""
# to allow code to work with Python 2 and 3


import os.path
import logging as log


class Demultiplex:
	"""
	It creates the command that will run the demultiplex
	steps. This step remove short sequences (i.e. adapters, polyA)
	from sequences
	"""

	def __init__(self, args):
		self.check_args(args)
		self._create_cmd()


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		self.wd = os.getcwd()
		self.params = args['params']
		self.execution = 1
		self.cmd = []
		self.cmd_file = ''
		if 'tmp_prefix' in args:
			self.tmp_prefix = args['tmp_prefix']
		else:
			if 'iter' in args:
				self.iter = args['iter']
				if args['iter'] == 'library':
					self.library = args['library']
					self.cmd_file = self.wd + '/' + __name__ + '_' + self.library + '_cmd.txt'
					self.tmp_prefix = self.library
				elif args['iter'] == 'sample':
					self.sample = args['sample']
					self.cmd_file = self.wd + '/' + __name__ + '_' + self.sample + '_cmd.txt'
					self.tmp_prefix = self.sample
				else:
					self.tmp_prefix = 'prefix'
			else:
				self.tmp_prefix = 'prefix'

		if 'middle' in args:
			self.middle = args['middle']
		else:
			self.middle = ''

		if 'min_qual' in args:
			self.min_qual = args['min_qual']

		if 'min_len' in args:
			self.min_len = args['min_len']

		if 'polyA' in args:
			self.polyA = args['polyA']

		if 'mid' in args:
			self.mid = args['mid']
		else:
			self.mid = ''

		if 'common' in args:
			self.common = args['common']
		else:
			self.common = ''

		if 'clean' in args:
			self.clean = bool(args['clean'])
		else:
			self.clean = False

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

		if 'adapters' in args:
			self.adapters = self._check_file(args['adapters'])
		else:
			self.adapters = ''

		if 'o1' in args:
			self.o1 = args['o1']

		if 'o2' in args:
			self.o2 = args['o2']

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
		"""
		Create command
		"""
		cmd = "#!/bin/bash\n"
		cmd += 'demultiplex.pl'
		if self.polyA:
			cmd += ' -polyA'
		if self.mid != '':
			for medium in self.mid:
				cmd += ' -index ' + medium + '=' + self.mid[medium]
		if self.common != '':
			cmd += ' -c ' + 'COMMON=' + self.common
		cmd += ' -1 ' + self.i1
		cmd += ' -2 ' + self.i2
		cmd += ' -tmp_prefix ' + self.tmp_prefix
		if self.middle != '':
			cmd += ' -middle ' + str(self.middle)
		cmd += ' -q ' + str(self.min_qual)
		cmd += ' -l ' + str(self.min_len)
		if self.clean:
			cmd += ' -clean '
		cmd += ' -a ' + self.adapters
		log.debug(cmd)
		self.cmd.append(cmd)


	def _check_file(self, file_path):
		"""
		Check that file exists
		"""
		try:
			open(file_path)
			return file_path
		except IOError:
			print('File not found ' + file_path)
