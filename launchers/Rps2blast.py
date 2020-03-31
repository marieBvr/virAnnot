# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import random
import string
from Bio import SeqIO

class Rps2blast:
	"""
	This module is part of virAnnot module
	This will select viral sequences from RPS results and
	run Blast on these sequences only.
	"""
	def __init__(self, args):
		self.execution = 1
		self.check_args(args)
		bcontigfile = os.getcwd() + '/' +str(args['sample']) + '/' + args['bcontigs']
		if not os.path.exists(bcontigfile):
			self.csv_to_fasta()
		self.cmd = []
		self.ssh_cmd = []
		if self.execution == 1:
			self.create_cmd()


	def create_cmd(self):
		"""
		Create command
		"""
		ssh_cmd = self._get_env_script()
		if self.server == 'genouest':
			exec_cmd = self._get_exec_script()
			cluster_cmd = self._get_clust_script()
		if self.server != 'enki':
			cmd = '#!/bin/bash\n'
			cmd += 'scp ' + self.bcontigs + ' ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
			cmd += ':' + self.params['servers'][self.server]['scratch']
			log.debug(cmd)
			self.cmd.append(cmd)
			fw = open(self.remote_cmd_file, mode='w')
			fw.write(ssh_cmd)
			fw.close()
			if self.server == 'genouest':
				fw = open(self.cluster_cmd_file, mode='w')
				fw.write(cluster_cmd)
				fw.close()
				fw = open(self.cluster_exec_cmd_file, mode='w')
				fw.write(exec_cmd)
				fw.close()
				cmd = 'ssh ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
				cmd += ' \'bash -s\' < ' + self.remote_cmd_file
				log.debug(cmd)
				self.cmd.append(cmd)
				cmd = 'scp ' + self.cluster_cmd_file + ' ' + self.params['servers'][self.server]['username']  + '@'
				cmd += self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir
				log.debug(cmd)
				self.cmd.append(cmd)
				cmd = 'ssh ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
				cmd += ' \'bash -s\' < ' + self.cluster_exec_cmd_file
				log.debug(cmd)
				self.cmd.append(cmd)
			else:
				cmd = 'ssh ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
				cmd += ' \'bash -s\' < ' + self.remote_cmd_file
				log.debug(cmd)
				self.cmd.append(cmd)
			cmd = 'scp ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress'] + ':'
			cmd += self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.out) + ' ' + self.wd
			log.debug(cmd)
			self.cmd.append(cmd)
		elif self.server == 'enki':
			self.cmd.append(ssh_cmd)


	def _get_exec_script(self):
		"""
		Script executed only on Genouest cluster
		"""
		exec_cmd = ''
		# exec_cmd += 'export SGE_ROOT="/usr/local/sge"\n'
		exec_cmd += '. /etc/profile.d/sge.sh\n'
		exec_cmd += 'echo "./' + os.path.basename(self.cluster_cmd_file)
		exec_cmd += '" | qsub -sync yes -V -wd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + ' -N ' + self.sample
		return exec_cmd


	def _get_clust_script(self):
		"""
		Script executed only on Genouest cluster
		"""
		clust_cmd = '. /softs/local/env/envpython-3.6.3.sh \n'
		clust_cmd += '. /softs/local/env/envblast-2.6.0.sh \n'
		clust_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + '\n'
		clust_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu ' + self.n_cpu + ' --tc ' + self.tc
		clust_cmd += ' -d ' + self.params['servers'][self.server]['db'][self.db]
		clust_cmd += ' -s ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.bcontigs)
		clust_cmd += ' --prefix ' + self.out_dir + ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + ' -r ' + ' --outfmt 5'
		clust_cmd += ' --max_target_seqs ' + self.max_target_seqs
		return clust_cmd


	def _get_env_script(self):
		"""
		This file must be in the cluster environment
		Depending on the cluster used, it load the
		correct environment and run the job command
		"""
		ssh_cmd = ''
		if self.server != 'enki':
			if self.server == 'genouest':
				ssh_cmd += '#!/bin/sh\n'
			else:
				ssh_cmd += '#!/bin/bash\n'
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
			if self.server == 'avakas':
				ssh_cmd += 'source ~/.bashrc' + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + "\n"
			ssh_cmd += 'mkdir ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"
			ssh_cmd += 'mv ' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.bcontigs) + ' ' + self.out_dir + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"
		if self.server != 'genouest':
			if self.server == 'genotoul':
				ssh_cmd += 'echo "'
			if self.server == "genologin":
				ssh_cmd += 'sbatch --mem=2G '
			ssh_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu ' + self.n_cpu
			ssh_cmd += ' --tc ' + self.tc + ' -d ' + self.params['servers'][self.server]['db'][self.db]
			if self.server != 'enki':
				ssh_cmd += ' -s ' + os.path.basename(self.bcontigs)
			else:
				ssh_cmd += ' -s ' + self.bcontigs

			ssh_cmd += ' --prefix ' + self.out_dir
			ssh_cmd += ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + ' -r ' + ' --outfmt 5'
			ssh_cmd += ' --max_target_seqs ' + self.max_target_seqs
			if self.server == 'genotoul':
				ssh_cmd += '"'
				ssh_cmd += ' | qsub -sync yes -V -wd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + ' -N ' + self.sample
		return ssh_cmd


	def csv_to_fasta(self):
		"""
		From rps csv results
		extract query id and generate fasta file
		"""
		fasta_file = self.icontigs # Input fasta file
		wanted_file = self.i # Input interesting sequence IDs
		result_file = self.bcontigs # Output fasta file

		wanted = set()
		with open(wanted_file) as f:
			for line in f:
				query_id = line.strip().split("\t")[0]
				query_length = line.strip().split("\t")[1]
				if query_length != "no_hit":
					wanted.add(query_id.replace("\"", "") )

		fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
		with open(result_file, "w") as f:
			for seq in fasta_sequences:
				if seq.id in wanted:
					SeqIO.write([seq], f, "fasta")


	def check_args(self, args=dict):
		"""
		Check if arguments are valid
		"""
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		accepted_type = ['tblastx', 'blastx', 'blastn', 'blastp']
		if 'i' in args:
			if os.path.exists(self.wd + '/' + args['i']):
				self.i = self.wd + '/' + args['i']
				self.execution = 1
			else:
				self.execution = 0
				log.critical('Input file do not exists.')
		if 'icontigs' in args:
			self.icontigs = self.wd + '/' + args['icontigs']
		else:
			self.icontigs = self.wd + '/' + self.sample + "_idba.scaffold.fa"
		if 'bcontigs' in args:
			self.bcontigs = self.wd + '/' + args['bcontigs']
		else:
			self.bcontigs = self.wd + '/' + self.sample + "_idba.scaffold.rps2bltx.fa"
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
		self.cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_rps2b_cmd.txt'
		self.remote_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_remote_rps2b_cmd.txt'
		self.cluster_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_cluster_rps2b_cmd.sh'
		self.cluster_exec_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_exec_rps2b_cmd.sh'
		self.random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
		self.out_dir = self.random_string + '_' + self.sample + '_' + self.type
