import argparse
import logging as log
import re
import os
import subprocess
import time
import glob
import shutil
import random
import warnings
from Bio import SeqIO

def main():
    args = _set_options()
    log_format = '%(asctime)s %(lineno)s %(levelname)-8s %(message)s'
    if(args.verbosity == 1):
        log.basicConfig(level=log.INFO,format=log_format)
    elif(args.verbosity == 3):
        log.basicConfig(level=log.DEBUG,format=log_format)
    wd = os.getcwd()
    out_dir = wd + '/' + args.prefix + '_split'
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    num_files = _split_fasta(args.seq,args.chunk,out_dir,args.random)
    jobid = _launch_jobs(num_files,out_dir,args.prefix,args.cluster,args.prog,args.db,args.n_cpu,args.tc,args.outfmt,args.max_target_seqs)
    _wait_job(args.cluster,jobid)
    if args.outfmt == 5:
        _concat_xml(out_dir,args.out)
    elif args.outfmt == 6:
        _concat_m8(out_dir,args.out)

    if args.clean:
        shutil.rmtree(out_dir)


def _concat_m8(path, out_file):
    files = glob.glob(path + '/' + '*.xml')
    for f in files:
        cmd = 'cat ' + f + ' >> ' + out_file
        os.system(cmd)


def _concat_xml(path, out_file):
    files = glob.glob(path + '/' + '*.xml')
    if len(files) == 1:
        shutil.move(files[0], out_file)
        return 0
    f_out = open(out_file, mode='w')
    h = None
    for f in files:
        h = open(f)
        header = h.readline()
        if not header:
            h.close()
            continue;
            warnings.warn("BLAST XML file %s was empty" % f)
        if header.strip() != '<?xml version="1.0"?>':
            f_out.write(header) #for diagnosis
            f_out.close()
            h.close()
            raise ValueError("%s is not an XML file!" % f)
        line = h.readline()
        header += line
        if line.strip() not in ['<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
                                '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">']:
            f_out.write(header) #for diagnosis
            f_out.close()
            h.close()
            raise ValueError("%s is not a BLAST XML file!" % f)
        while True:
            line = h.readline()
            if not line:
                f_out.write(header) #for diagnosis
                f_out.close()
                h.close()
                log.critical('BLAST XML file %s ended prematurely' % f)
                raise ValueError("BLAST XML file %s ended prematurely" % f)
            header += line
            if "<Iteration>" in line:
                break
            if len(header) > 10000:
                #Something has gone wrong, don't load too much into memory!
                #Write what we have to the merged file for diagnostics
                f_out.write(header)
                f_out.close()
                h.close()
                raise ValueError("BLAST XML file %s has too long a header!" % f)
        if "<BlastOutput>" not in header:
            f_out.close()
            h.close()
            log.critical("%s is not a BLAST XML file:\n%s\n..." % (f, header))
            raise ValueError("%s is not a BLAST XML file:\n%s\n..." % (f, header))
        if f == files[0]:
            f_out.write(header)
        else:
            f_out.write("    <Iteration>\n")
        for line in h:
            if "</BlastOutput_iterations>" in line:
                break
            #TODO - Increment <Iteration_iter-num> and if required automatic query names
            #like <Iteration_query-ID>Query_3</Iteration_query-ID> to be increasing?
            f_out.write(line)
        h.close()
    f_out.write("  </BlastOutput_iterations>\n")
    f_out.write("</BlastOutput>\n")
    f_out.close()


def _wait_job(cluster=str, jobid=int):
    qstat_cmd = _get_qstat_cmd(cluster,jobid)
    log.info('Waiting job array ' + jobid)
    pipes = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipes.communicate()
    stdoutdata = stdout.decode("utf-8").split("\n")
    while len(stdoutdata) > 2:
        log.debug('job ' + jobid + ' still runing...')
        time.sleep(20)
        pipes = subprocess.Popen(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipes.communicate()
        stdoutdata = stdout.decode("utf-8").split("\n")
        log.debug(stdoutdata)


def _get_qstat_cmd(cluster=str, jobid=int) :
    qstat_cmd = ''
    if cluster == 'enki':
        qstat_cmd = ['qstat', '-j',str(jobid)]
    elif cluster == 'genologin':
        qstat_cmd = ['squeue', '-j',str(jobid)]
    elif cluster == 'genouest':
        qstat_cmd = ['qstat', '-j',str(jobid)]
    elif cluster == 'curta':
        qstat_cmd = ['qstat', '-i',str(jobid)]
    else:
        log.critical('unknown cluster.')
    return qstat_cmd


def _launch_jobs(n_f, o_d, prefix, cluster, prog, db, n_cpu, tc, outfmt, max_target_seqs):
    script = _load_script(cluster,prog)
    blt_script = o_d + '/' + 'blast_script.sh'
    fw = open(blt_script, mode='w')
    fw.write(script % (prog, o_d, db, o_d, 0.001, outfmt, max_target_seqs, n_cpu))
    fw.close()
    qsub_cmd, job_regex = _get_qsub_cmd(cluster,n_f,o_d,n_cpu,tc,blt_script)
    log.info(qsub_cmd)
    pipes = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipes.communicate()
    stdoutdata = stdout.decode("utf-8")
    log.debug(stdoutdata)
    # stdoutdata = subprocess.getoutput(qsub_cmd)
    # log.debug(stdoutdata)
    p = re.compile(job_regex)
    m = p.match(stdoutdata)
    try:
        jobid = m.group(1)
    except AttributeError:
        jobid = m.group()
    log.debug('job launch with jobid ' + jobid + '.')
    time.sleep(10)
    return jobid

def _get_qsub_cmd(cluster=str, n_f=str, o_d=str, n_cpu=int, tc=int, blt_script=str):
    qsub_cmd = ''
    job_regex = ''
    if cluster == 'enki':
        # qsub_cmd = 'qsub -V -t 1-' + str(n_f+1) + ' -wd ' + o_d + ' -pe multithread ' + str(n_cpu) + ' ' + blt_script
        qsub_cmd = ['qsub', '-V', '-t', '1-' + str(n_f+1), '-tc', str(tc), '-wd', o_d, '-pe', 'multithread', str(n_cpu), blt_script]
        job_regex = '^Your job-array (\d+)\.1-\d+'
    elif cluster == 'curta':
        qsub_cmd = ['sbatch', '--export=ALL', '--array=1-' + str(n_f+1), '-D' , o_d, '--time=250:00:00', '--mem=14G', '--nodes=1 --ntasks=',
            str(n_cpu), blt_script]
        # qsub_cmd = 'qsub -V -t 1-' + str(n_f+1) + ' -d ' + o_d + ' -l walltime=18:00:00 -l nodes=1:ppn=' + str(n_cpu) + ' ' + blt_script
        job_regex = '^Submitted batch job (\d+)'
    elif cluster == 'genouest':
        qsub_cmd = ['qsub', '-V', '-t', '1-' + str(n_f+1), '-tc', str(tc), '-wd', o_d,'-l', 'mem=8G', '-l', 'h_vmem=12G', '-pe', 'parallel_smp',
            str(n_cpu), blt_script]
        job_regex = '^Your job-array (\d+)\.1-\d+:\d+ \(\"blast_script\.sh\"\) has been submitted'
    elif cluster == 'genologin':
        qsub_cmd = ['sbatch','--export=ALL', '--array=1-' + str(n_f+1), '--ntasks-per-node=' + str(tc), '-D', o_d, '--mem=20G',
            '--ntasks=' + str(n_cpu), blt_script]
        job_regex = '^Submitted batch job (\d+)'
        # job_regex = '^Waiting job array (\d+)'

    else:
        log.critical('unknown cluster.')
    log.debug(qsub_cmd)
    return qsub_cmd, job_regex


def _split_fasta(fasta,chunk,directory,rand):
    record_iter = SeqIO.parse(open(fasta),'fasta')
    if rand:
        record_iter = sorted(record_iter, key=lambda k : random.random())
    for i, batch in enumerate(batch_iterator(iter(record_iter), chunk)):
        filename = directory + '/' + "group_%i.fa" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
    log.info(fasta + ' file splited in ' + str(i+1) + ' part in ' + directory + '.')
    return i


def batch_iterator(iterator, batch_size):
    entry = True
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break
            batch.append(entry)
        if batch:
            yield batch


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--seq',help='The fasta sequence file.',action='store',type=str,required=True)
    parser.add_argument('-c','--cluster',help='The cluster name.',action='store',type=str,required=True,default='avakas',choices=['enki','genologin','genouest', 'curta'])
    parser.add_argument('-n','--num_chunk',dest='chunk',help='The number of chunks to split the fasta in.',action='store',type=int,default=100)
    parser.add_argument('--tc',dest='tc',help='The number of concurent jobs to launch on SGE servers.',action='store',type=int,default=100)
    parser.add_argument('--n_cpu',help='The number of cpu cores to use per job.',action='store',type=int,default=5)
    parser.add_argument('-p','--prog',help='The Blast program to use.',action='store',type=str,default='blastx',choices=['blastx','blastn','blastp','tblastx','rpstblastn'])
    parser.add_argument('-d','--db',help='The Blast database to use.',action='store',type=str,default='nr')
    parser.add_argument('-o','--out',help='Output XML file.',action='store',type=str,default='blast-out.xml')
    parser.add_argument('--outfmt',help='Output Blast format. 5: XML, 6: m8',action='store',type=int,default=5,choices=[5,6])
    parser.add_argument('--max_target_seqs',help='Maximum number of aligned sequences to keep.',action='store',type=int,default=5)
    parser.add_argument('--prefix',help='Directory prefix to store splited files.',action='store',type=str,default='split')
    parser.add_argument('--clean',help='Delete blast directory.',action='store_true')
    parser.add_argument('-r','--random',help='Randomly split sequences to balance load between jobs.',action='store_true')
    parser.add_argument('-v','--verbosity',help='Verbose level', action='store',type=int,choices=[1,2,3,4],default=1)
    args = parser.parse_args()
    return args


def _load_script(cluster, prog):
    script_sge = ''
    if cluster == 'enki':
        script_sge = "#!/bin/sh\n%s -query %s/group_$SGE_TASK_ID.fa -db %s -out %s/group_$SGE_TASK_ID.xml -evalue %f -outfmt %d -max_target_seqs %d -parse_deflines -num_threads %i"
    elif cluster == 'genouest':
        script_sge = "#!/bin/sh\n"
        script_sge += ".  /softs/local/env/envblast-2.6.0.sh\n"
        script_sge += ". /softs/local/env/envpython-3.6.3.sh\n"
        script_sge += "%s -query %s/group_$SGE_TASK_ID.fa -db %s -out %s/group_$SGE_TASK_ID.xml -evalue %f -outfmt %d -max_target_seqs %d -parse_deflines -num_threads %i"
    elif cluster == 'genologin':
        script_sge = "#!/bin/sh\n%s -query %s/group_$SLURM_ARRAY_TASK_ID.fa -db %s -out %s/group_$SLURM_ARRAY_TASK_ID.xml -evalue %f -outfmt %d -max_target_seqs %d -parse_deflines -num_threads %i"
    elif  cluster == 'curta':
        script_sge = "#!/bin/bash\n%s -query %s/group_$SLURM_ARRAY_TASK_ID.fa -db %s -out %s/group_$SLURM_ARRAY_TASK_ID.xml -evalue %f -outfmt %d -max_target_seqs %d -parse_deflines -num_threads %i"
    return script_sge


if __name__ == "__main__":
    main()
