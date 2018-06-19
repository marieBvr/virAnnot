import os.path
from subprocess import call
import logging as log
import random
import string
import re

"""
title           :Diamod2blast.py
description     :This will select viral sequences from DIAMOND results and 
				  run Blast on these sequences only.
author          :Marie Lefebvre
date            :june-07-2018
version         :0.1
usage           :python3.4 Diamod2blast.py
notes           :
python_version  :3.4  
==============================================================================
"""


class Diamond2blast:
	def __init__ (self, args):
		self.check_args(args)
		self.csv_to_fasta()
		self.cmd = []
		self.ssh_cmd = []
		if self.execution == 1:
			self.create_cmd()



	def create_cmd (self):
		ssh_cmd = self._get_exec_script()
		if self.server != 'enki' :
			if self.server == 'avakas':
				cmd = 'scp ' + self.contigs + ' ' + self.username + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch']
			else:
				cmd = 'scp ' + self.contigs + ' ' + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch']
			log.debug(cmd)
			self.cmd.append(cmd)
			fw =  open(self.remote_cmd_file, mode='w')
			fw.write(ssh_cmd)
			fw.close()
			if self.server != 'avakas':
				cmd = 'ssh ' + self.username + '@' + self.params['servers'][self.server]['adress'] + ' \'bash -s\' < ' + self.remote_cmd_file
				log.debug(cmd)
				self.cmd.append(cmd)
				cmd = 'scp ' + self.username + '@' + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.out) + ' ' + self.wd
				log.debug(cmd)
				self.cmd.append(cmd)
			else:
				cmd = 'ssh ' + self.username + '@' + self.params['servers'][self.server]['adress'] + ' \'bash -s\' < ' + self.remote_cmd_file
				log.debug(cmd)
				self.cmd.append(cmd)
				cmd = 'scp ' + self.username + '@' + self.params['servers'][self.server]['adress'] + ':' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.out) + ' ' + self.wd
				log.debug(cmd)
				self.cmd.append(cmd)
		elif self.server == 'enki':
			self.cmd.append(ssh_cmd)


	def _get_exec_script (self):
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
			if self.server == 'avakas':
				ssh_cmd += 'source ~/.bashrc' + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + "\n"
			ssh_cmd += 'mkdir ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"
			ssh_cmd += 'mv ' + self.params['servers'][self.server]['scratch'] + '/' + os.path.basename(self.contigs) + ' ' + self.out_dir + "\n"
			ssh_cmd += 'cd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + "\n"

		if self.server == 'genotoul':
			ssh_cmd += 'echo "'
		# elif self.server == "genologin":
		# 	ssh_cmd += 'echo "'
		ssh_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu ' + self.n_cpu + ' --tc ' + self.tc + ' -d ' + self.params['servers'][self.server]['db'][self.db]
		if self.server != 'enki':
			ssh_cmd += ' -s ' + os.path.basename(self.contigs)
		else:
			ssh_cmd += ' -s ' + self.contigs

		ssh_cmd += ' --prefix ' + self.out_dir
		ssh_cmd += ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + ' -r ' + ' --outfmt 5'
		ssh_cmd += ' --max_target_seqs ' + self.max_target_seqs
		if self.server == 'genotoul':
			ssh_cmd += '"'
			ssh_cmd += ' | qsub -sync yes -V -wd ' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + ' -N ' + self.sample
		# elif self.server == 'genologin':
		# 	ssh_cmd += '"'
		# 	ssh_cmd += ' | sbatch -W --export=ALL --workdir=' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + ' -J ' + self.sample
		return ssh_cmd


	def csv_to_fasta (self):
		#File input
		fileInput = open(self.i, "r")

		#File output
		fileOutput = open(self.contigs, "w")

		#Seq count
		count = 1 ;

		#Loop through each line in the input file
		log.debug("Converting to FASTA...")
		for strLine in fileInput:
		    reg = re.compile(r"^\"Virus.*")
		    taxo = strLine.split("\t")[14]
		    #Test if viral
		    if re.match(reg, taxo) != None:
		    	seq = strLine.split("\t")[15]
		    	#Write sequence header
		    	fileOutput.write(">" + self.sample + "_d2b_" + str(count) + "\n")
		    	#Write sequence
		    	fileOutput.write(seq.replace('\"',''))
		    	count = count + 1
		log.debug("FASTA conversion done for" + self.sample)

		#Close the input and output file
		fileInput.close()
		fileOutput.close()


	def check_args (self, args: dict):
		if 'sample' in args:
			self.sample = str(args['sample'])
		self.wd = os.getcwd() + '/' + self.sample
		accepted_type = ['tblastx','blastx','blastn','blastp','rpstblastn']
		if 'i' in args:
			if os.path.exists(self.wd + '/' + args['i']):
				self.i = self.wd + '/' + args['i']
				self.execution = 1;
			else:
				self.execution = 0;
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
				log.critical(arg['db'] + ' not defined in parameters file')
			else:
				self.db = args['db']
		else:
			log.critical('You must provide a database name.')
		self.cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_d2b_cmd.txt'
		self.remote_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_remote_d2b_cmd.txt'
		self.random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
		self.out_dir = self.random_string + '_' + self.sample + '_' + self.type
