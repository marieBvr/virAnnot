"""
    This class is a part of the virAnnot module
    ============
    Authors: Sebastien Theil, Marie Lefebvre
"""

# to allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log
import random
import string

class Blast:

    def __init__(self, args):
        self.check_args(args)
        self.cmd = []
        self.ssh_cmd = []
        if self.execution == 1:
            self.create_cmd()


    def create_cmd(self):
        #
        # Send contig files and scripts to cluster to execute command
        #
        ssh_cmd = self.get_exec_script()
        if self.server != 'enki':
            cmd = '#!/bin/bash\n'
            cmd += 'scp ' + self.contigs + ' ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
            cmd += ':' + self.params['servers'][self.server]['scratch']
            log.debug(cmd)
            self.cmd.append(cmd)
            if self.server == 'genouest':
                cmd = 'scp ' + self.genouest_cmd_file + ' ' + self.params['servers'][self.server]['username']
                cmd += '@' + self.params['servers'][self.server]['adress']
                cmd += ':' + self.params['servers'][self.server]['scratch']
                log.debug(cmd)
                self.cmd.append(cmd)
            fw = open(self.remote_cmd_file, mode='w')
            fw.write(ssh_cmd)
            fw.close()
            if self.server != 'curta':
                cmd = 'ssh ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
                cmd += ' \'bash -s\' < ' + self.remote_cmd_file
                log.debug(cmd)
                self.cmd.append(cmd)
                cmd = 'scp ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress'] + ':'
                cmd += self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.out) + ' ' + self.wd
                log.debug(cmd)
                self.cmd.append(cmd)
            else:
                cmd = 'ssh ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
                cmd += ' \'bash -s\' < ' + self.remote_cmd_file
                log.debug(cmd)
                self.cmd.append(cmd)
                cmd = 'scp ' + self.params['servers'][self.server]['username'] + '@' + self.params['servers'][self.server]['adress']
                cmd += ':' + self.params['servers'][self.server]['scratch'] + '/' + self.out_dir + '/' + os.path.basename(self.out) + ' ' + self.wd
                log.debug(cmd)
                self.cmd.append(cmd)
        elif self.server == 'enki':
            self.cmd.append(ssh_cmd)


    def get_exec_script(self):
        #
        # Set cluster environment and launch blast_launch.py
        # For genouest cluster, need an additional file
        #
        ssh_cmd = ''
        log.info("get exec")
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
            if self.server == "enki":
                ssh_cmd += "#!/bin/sh\n"
            ssh_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu 8 --tc ' + self.tc
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
        #
        # Create cmd file to launch blast_launch.py on genouest cluster
        #
        go_cmd = '#!/bin/bash' + "\n"
        go_cmd += '. /local/env/envconda.sh' + "\n"
        go_cmd += 'conda activate ~/blast_env' + "\n"
        go_cmd += 'blast_launch.py -c ' + self.server + ' -n ' + self.num_chunk + ' --n_cpu 8 --tc ' + self.tc
        go_cmd += ' -d ' + self.params['servers'][self.server]['db'][self.db] + ' -s ' + os.path.basename(self.contigs)
        go_cmd += ' --prefix ' + self.out_dir
        go_cmd += ' -p ' + self.type + ' -o ' + os.path.basename(self.out) + ' -r ' + ' --outfmt 5'
        go_cmd += ' --max_target_seqs ' + self.max_target_seqs
        fw = open(self.genouest_cmd_file, mode='w')
        fw.write(go_cmd)
        fw.close()


    def check_args(self, args=dict):
        if 'sample' in args:
            self.sample = str(args['sample'])
        self.wd = os.getcwd() + '/' + self.sample
        accepted_type = ['tblastx', 'blastx', 'blastn', 'blastp', 'rpstblastn']
        if 'contigs' in args:
            if os.path.exists(self.wd + '/' + args['contigs']):
                self.contigs = self.wd + '/' + args['contigs']
                self.execution = 1
            else:
                self.execution = 0
                log.critical('Input fasta file do not exists. ' + self.wd + '/' + args['contigs'])
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
        if 'db' in args:
            if args['db'] not in self.params['servers'][self.server]['db']:
                log.critical(args['db'] + ' not defined in parameters file')
            else:
                self.db = args['db']
        else:
            log.critical('You must provide a database name.')
        self.cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_blast_cmd.txt'
        self.remote_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_remote_blast_cmd.txt'
        self.genouest_cmd_file = self.wd + '/' + self.sample + '_' + self.type + '_' + self.db + '_genouest_cmd.txt'
        self.random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
        self.out_dir = self.random_string + '_' + self.sample + '_' + self.type
