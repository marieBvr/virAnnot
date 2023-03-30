"""
    This class is a part of the virAnnot module
    ============
    Authors: Sebastien Theil, Marie Lefebvre
"""

# To allow code to work with Python 2 and 3
from __future__ import print_function   # print is a function in python3
from __future__ import unicode_literals # avoid adding "u" to each string
from __future__ import division # avoid writing float(x) when dividing by x

import os.path
import logging as log

class Rps2ecsv:
    """
    Convert RPSblast results (XML) to CSV file
    """

    def __init__(self, args):
        self.check_args(args)
        self.cmd = []
        self.create_cmd()


    def create_cmd(self):
        cmd = '#!/bin/sh\n'
        cmd += 'rps2ecsv.py'
        cmd += ' -b ' + self.b
        cmd += ' -e ' + str(self.evalue)
        cmd += ' -o ' + self.out
        self.cmd.append(cmd)


    def check_args(self, args=dict):
        self.execution = 1
        if 'sample' in args:
            self.sample = str(args['sample'])
        self.wd = os.getcwd() + '/' + self.sample
        self.cmd_file = self.wd + '/' + self.sample + '_rps2ecsv_cmd.txt'
        if 'contigs' in args:
            self.contigs = self.wd + '/' + args['contigs']
        if 'sge' in args:
            self.sge = bool(args['sge'])
        else:
            self.sge = False
        if 'out' in args:
            self.out = args['out']
        if 'params' in args:
            self.params = args['params']
        if 'evalue' in args:
            self.evalue = args['evalue']
        else:
            self.evalue = 0.0001
        if 'b' in args:
            if os.path.exists(self.wd + '/' + args['b']):
                self.b = self.wd + '/' + args['b']
            else:
                self.b = ''
                self.execution = 0
        else:
            log.critical('You must provide a blast rps file.')
        if 'n_cpu' in args:
            self.n_cpu = str(args['n_cpu'])
        else:
            self.n_cpu = '1'
