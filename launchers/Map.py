# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log

class Map:
	"""
	This class is a part of the virAnnot module
	"""

	def __init__(self, args):
		self.check_args(args)
		self.cmd = []
		self._create_cmd()

	def _create_cmd(self):
		"""
		Create command
		"""
		cmd = '#!/bin/sh\n bowtie2-build ' + self.contigs + ' ' + self.contigs
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = 'bowtie2 -p ' + self.n_cpu + ' -x ' + self.contigs + ' -1 ' + self.i1 + ' -2 ' + self.i2
		if self.ising != '':
			cmd += ' -U ' + self.ising
		# cmd = 'fa2fq ' + self.ising + ' ' + os.path.splitext(self.ising)[0] + '.fastq'
		# log.debug(cmd)
		# self.cmd.append(cmd)
		# cmd = 'bowtie2 -p ' + self.n_cpu + ' -x ' + self.contigs + ' -U ' + os.path.splitext(self.ising)[0] + '.fastq'
		cmd += ' | ' + self.params['bin']['samtools'] + ' view -bS - > ' + self.bam
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = self.params['bin']['samtools'] + ' sort ' + self.bam + ' -o ' + os.path.splitext(self.bam)[0] + '.sort.bam'
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = self.params['bin']['samtools'] + ' index ' + os.path.splitext(self.bam)[0] + '.sort.bam'
		log.debug(cmd)
		self.cmd.append(cmd)
		cmd = 'readPerContig.pl -r ' + self.contigs + ' -o ' + self.rn + ' -p bowtie -b ' + os.path.splitext(self.bam)[0] + '.sort.bam'
		log.debug(cmd)
		self.cmd.append(cmd)


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		self.cmd_file = self.wd + '/' + self.sample + '_map_cmd.txt'
		self.execution = 1
		if 'contigs' in args:
			self.contigs = self.wd + '/' + args['contigs']
		if 'i1' in args:
			self.i1 = self.check_seq_format(self.wd + '/' + args['i1'])
		else:
			log.critical('Need r1 file.')
			self.execution = 0
		if 'i2' in args:
			self.i2 = self.check_seq_format(self.wd + '/' + args['i2'])
		else:
			log.critical('Need r2 file.')
			self.execution = 0
		if 'ising' in args:
			self.ising = self._check_file(self.wd + '/' + args['ising'])
		else:
			self.ising = ''
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			log.debug('n_cpu option not found. default 1')
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'bam' in args:
			self.bam = self.wd + '/' + args['bam']
		if 'rn' in args:
			self.rn = self.wd + '/' + args['rn']
		if 'params' in args:
			self.params = args['params']


	def launch(self):
		if self.sge:
			fw = open(self.cmd_file, mode='w')
			for el in self.cmd:
				fw.write(el + "\n")
			fw.close()
			qsub_call = "qsub -wd " + self.wd + " -V -N " + self.sample + '_map' + ' -pe multithread ' + self.n_cpu + ' ' + self.cmd_file
			log.debug(qsub_call)
			os.system(qsub_call)
		else:
			for el in self.cmd:
				os.system(el)


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

	def check_seq_format(self, in_file):
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
