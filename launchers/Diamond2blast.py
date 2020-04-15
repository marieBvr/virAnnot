"""
title           :Diamod2blast.py
author          :Marie Lefebvre
date            :june-07-2018
version         :0.1
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import random
import string
import re
from Blast import Blast


class Diamond2blast:
	"""
	This module is part of virAnnot module
	This will select viral sequences from DIAMOND results and
	run Blast on these sequences only.
	"""
	def __init__(self, args):
		self.check_args(args)
		contigfile = os.getcwd() + '/' +str(args['sample']) + '/' + args['contigs']
		if not os.path.exists(contigfile):
			self.csv_to_fasta()
		self.cmd = []
		self.ssh_cmd = []
		if self.execution == 1:
			Blast.create_cmd(self)


	def get_exec_script(self):
		#
		# Set cluster environment and launch blast_launch.py
		# For genouest cluster, need an additional file
		#
		ssh_cmd = ''
		if self.server != 'enki':
			ssh_cmd += 'if [ -f ~/.bashrc ]; then' + "\n"
			ssh_cmd += 'source ~/.bashrc' + "\n"
			ssh_cmd += 'echo bashrc loaded' + "\n"
			ssh_cmd += 'elif [ -f ~/.profile ]; then' + "\n"
			ssh_cmd += 'source ~/.profile' + "\n"
			ssh_cmd += 'echo profile loaded' + "\n"
			ssh_cmd += 'elif [ -f ~/.bash_profile ]; then' + "\n"
			ssh_cmd += 'source /etc/profile' + "\n"
			ssh_cmd += 'source ~/.bash_profile' + "\n"
			ssh_cmd += 'echo bash_profile loaded' + "\n"
			ssh_cmd += 'else' + "\n"
			ssh_cmd += 'echo "source not found."' + "\n"
			ssh_cmd += 'fi' + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + "\n"
			ssh_cmd += 'mkdir ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"
			ssh_cmd += 'mv ' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.contigs) + ' ' + self.out_dir + "\n"
			if self.server == 'genouest':
				ssh_cmd += 'mv ' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.genouest_cmd_file) + ' ' + self.out_dir + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"
		if self.server == 'genouest':
			self.create_genouest_script()
			ssh_cmd += 'sbatch ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir
			ssh_cmd += '/' + os.path.basename(self.genouest_cmd_file)
		else:
			if self.server == "genologin":
				ssh_cmd += 'sbatch --mem=2G '
			ssh_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu ' + self.n_cpu + ' --tc ' + self.tc
			ssh_cmd += ' -d ' + self.params['servers'][self.server]['db'][self.db]
			if self.server != 'enki':
				ssh_cmd += ' -s ' + os.path.basename(self.contigs)
			else:
				ssh_cmd += ' -s ' + self.contigs
			ssh_cmd += ' --prefix ' + self.out_dir
			ssh_cmd += ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + ' -r ' + ' --outfmt 5'
			ssh_cmd += ' --max_target_seqs ' + self.max_target_seqs
		return ssh_cmd


	def create_genouest_script(self):
		Blast.create_genouest_script(self)


	def csv_to_fasta(self):
		"""
		From diamond csv results
		extract Viral sequences and
		generate fasta file
		"""
		#File input
		file_input = open(self.i, "r")
		#File output
		file_output = open(self.contigs, "w")
		#Seq count
		count = 1

		#Loop through each line in the input file
		log.debug("Converting to FASTA...")
		for strLine in file_input:
			reg = re.compile(r"^\"Virus.*")
			taxo = strLine.split("\t")[14]
			query_id = strLine.split("\t")[1]
			#Test if viral
			if re.match(reg, taxo) != None:
				seq = strLine.split("\t")[15]
				#Write sequence header
				file_output.write(">" + query_id.replace("\"", "") + "\n")
				#Write sequence
				file_output.write(seq.replace('\"', ''))
				count = count + 1
		log.debug("FASTA conversion done for " + self.sample)

		#Close the input and output file
		file_input.close()
		file_output.close()


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		accepted_type = ['tblastx', 'blastx', 'blastn', 'blastp', 'rpstblastn']
		if 'i' in args:
			if os.path.exists(self.wd + '/' + args['i']):
				self.i = self.wd + '/' + args['i']
				self.execution = 1
			else:
				self.execution = 0
				log.critical('Input file do not exists.')
		if 'contigs' in args:
			self.contigs = self.wd + '/' + args['contigs']
		else:
			self.contigs = self.wd + '/' + self.sample + "_idba.scaffold.dmdx2bltx.fa"
		if 'type' in args:
			if args['type'] in accepted_type:
				self.type = args['type']
			else:
				log.critical('Wrong blast type. ' + accepted_type)
		else:
			log.critical('Blast type is mandatory.')
		if 'n_cpu' in args:
			self.n_cpu = str(args['n_cpu'])
		else:
			log.debug('n_cpu option not found. default 1')
			self.n_cpu = '1'
		if 'sge' in args:
			self.sge = bool(args['sge'])
		else:
			self.sge = False
		if 'tc' in args:
			self.tc = str(args['tc'])
		else:
			self.tc = '5'
		if 'max_target_seqs' in args:
			self.max_target_seqs = str(args['max_target_seqs'])
		else:
			self.max_target_seqs = '5'

		if 'num_chunk' in args:
			self.num_chunk = str(args['num_chunk'])
		else:
			self.num_chunk = '100'
		if 'out' in args:
			self.out = args['out']
		if 'params' in args:
			self.params = args['params']
		if 'server' in args:
			self.server = args['server']
		if 'username' in args['params']['servers'][self.server]:
			self.username = args['params']['servers'][self.server]['username']
		else:
			log.critical('No username defined for cluster.')
		if 'db' in args:
			if args['db'] not in self.params['servers'][self.server]['db']:
				log.critical(args['db'] + ' not defined in parameters file')
			else:
				self.db = args['db']
		else:
			log.critical('You must provide a database name.')

		self.cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_d2b_cmd.txt'
		self.remote_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_remote_d2b_cmd.txt'
		self.genouest_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_genouest_d2b_cmd.txt'
		self.random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
		self.out_dir = self.random_string + '_' + self.sample + '_' + self.type

