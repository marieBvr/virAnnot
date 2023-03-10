#!/usr/bin/env python3.6
###!/usr/bin/python3.6

import argparse
import logging as log
import csv
import importlib
import re
import sys
import os, shutil
import time
import yaml


"""
Main function
"""
def main():
	args = _set_options()
	_set_log_level(args.verbosity)

	params = _read_yaml_file(args.param)
	steps = _read_yaml_file(args.step)
	maps = _read_map_file(args.map)
	if args.name_step == 'init':
		log.info('Init directory and move files...')
		_create_folders(maps)
	elif(args.name_step in steps):
		log.info('Launching step ' + args.name_step)
		start_time = time.time()
		_launch_step(args.name_step,steps,maps,params)
		log.info("--- %s seconds ---" %(time.time() - start_time))
	else:
		log.critical('This step is not present in the step file.')

def _set_log_level(verbosity):
	if verbosity == 1:
		log_format = '%(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.INFO,format=log_format)
	elif verbosity == 3:
		log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.DEBUG,format=log_format)

def _create_folders(maps):
	wd = os.getcwd()
	for i in range(0,len(maps)):
		if 'file' in maps[i]:
			if not os.path.exists(maps[i]['file']):
				log.critical(maps[i]['file'] + ' not found.')
				sys.exit(1)
			else:
				if not os.path.exists(wd + '/' + maps[i]['SampleID']):
					os.mkdir(wd + '/' + maps[i]['SampleID'])
					shutil.move(wd + '/' + maps[i]['file'],wd + '/' + maps[i]['SampleID'])
		elif 'file1' in maps[i] and 'file2' in maps[i]:
			if not os.path.exists(maps[i]['file1']):
				log.critical(maps[i]['file1'] + ' not found.')
				sys.exit(1)
			if not os.path.exists(maps[i]['file2']):
				log.critical(maps[i]['file2'] + ' not found.')
				sys.exit(1)
			else:
				if not os.path.exists(wd + '/' + maps[i]['SampleID']):
					os.mkdir(wd + '/' + maps[i]['SampleID'])
					shutil.move(wd + '/' + maps[i]['file1'],wd + '/' + maps[i]['SampleID'])
					shutil.move(wd + '/' + maps[i]['file2'],wd + '/' + maps[i]['SampleID'])
				else:
					shutil.move(wd + '/' + maps[i]['file1'],wd + '/' + maps[i]['SampleID'])
					shutil.move(wd + '/' + maps[i]['file2'],wd + '/' + maps[i]['SampleID'])
		else:
			log.critical('None of the required file column found.')
			sys.exit(1)


def _launch_step(s_n,s,m,p):
	module_name = s_n.split('_')[0]
	if module_name == 'Getresults':
		_launch_getresults(s_n,s,m,p,module_name)
	else:
		if('iter' in s[s_n]):
			if(s[s_n]['iter'] == 'library'):
				 _launch_by_library(s_n,s,m,p,module_name)
			elif(s[s_n]['iter'] == 'sample'):
				_launch_by_sample_id(s_n,s,m,p,module_name)
			elif(s[s_n]['iter'] == 'global'):
				_launch_global(s_n,s,m,p,module_name)
			else:
				log.critical('iter options must be library, sample or global')
				sys.exit(1)
		else:
			_launch_by_sample_id(s_n,s,m,p,module_name)


def _launch_getresults(s_n,s,m,p,module_name):
	args={}
	args['global_dir']=[]
	args['sample_dir']={}
	args['global_files']=[]
	args['sample_files']={}
	for key in s[s_n]:
		if 'global_dir' in key:
			args['global_dir'].append(s[s_n][key])
		if 'global_file' in key:
			args['global_files'].append(s[s_n][key])
		if key == 'out':
			args['out'] = s[s_n][key]
		if 'sample' in key:
			for i in range(0,len(m)):
				tmp = _replace_sample_name(s[s_n],m[i])
				if 'dir' in key:
					if m[i]['SampleID'] not in args['sample_dir']:
						args['sample_dir'][m[i]['SampleID']]=[]
					args['sample_dir'][m[i]['SampleID']].append(tmp[key])
				if 'file' in key:
					if m[i]['SampleID'] not in args['sample_files']:
						args['sample_files'][m[i]['SampleID']]=[]
					args['sample_files'][m[i]['SampleID']].append(tmp[key])
	_launch_module(args,module_name)

def _launch_global(s_n,s,m,p,module_name):
	global_input = {}
	global_input['args'] = {}
	for i in range(0,len(m)):
		tmp = _replace_sample_name(s[s_n],m[i])
		global_input['args'][m[i]['SampleID']] = {}
		for j in tmp:
			if j == 'out':
				if 'out' not in global_input:
					global_input['out'] = tmp[j]
					continue
			elif j == 'outdir':
				if 'outdir' not in global_input:
					global_input['outdir'] = tmp[j]
					continue
			elif j == 'sge':
				if 'sge' not in global_input:
					global_input['sge'] = tmp[j]
					continue
			elif j == 'data':
				if  'data' not in global_input:
					global_input['data'] = tmp[j]
					continue
			elif j == 'r':
				if 'r' not in global_input:
					global_input['r'] = tmp[j]
					continue
			elif j == 'c':
				if 'c' not in global_input:
					global_input['c'] = tmp[j]
					continue
			elif j == 'viral_portion':
				if 'viral_portion' not in global_input:
					global_input['viral_portion'] = tmp[j]
					continue
			elif j == 'iter':
				if 'iter' not in global_input:
					global_input['iter'] = tmp[j]
					continue
			elif j == 'min_prot':
				if 'min_prot' not in global_input:
					global_input['min_prot'] = tmp[j]
					continue
			elif j == 'blast_db':
				if 'blast_db' not in global_input:
					global_input['blast_db'] = tmp[j]
					continue
			elif j == 'blast_type':
				if 'blast_type' not in global_input:
					global_input['blast_type'] = tmp[j]
					continue
			elif j == 'perc':
				if 'perc' not in global_input:
					global_input['perc'] = tmp[j]
					continue
			elif j == 'rps_folder':
				if 'rps_folder' not in global_input:
					global_input['rps_folder'] = tmp[j]
					continue
			else:
				global_input['args'][m[i]['SampleID']][j] = tmp[j]
	global_input['params'] = p
	_launch_module(global_input,module_name)


def _launch_by_sample_id(s_n,s,m,p,module_name):
	for i in range(0,len(m)):
		tmp = _replace_sample_name(s[s_n],m[i])
		tmp['sample'] = m[i]['SampleID']
		tmp['params'] = p
		_launch_module(tmp,module_name)


def _launch_by_library(s_n,s,m,p,module_name):
	args = {}
	for i in range(0,len(m)):
		tmp = _replace_sample_name(s[s_n],m[i])
		if(m[i]['library'] not in args):
			args[m[i]['library']] = {}
		if('file1' not in args[m[i]['library']] and 'i1' not in s[s_n]):
			args[m[i]['library']]['i1'] = m[i]['file1']
		else:
			args[m[i]['library']]['i1'] = tmp['i1']
		if('file2' not in args[m[i]['library']] and 'i2' not in s[s_n]):
			args[m[i]['library']]['i2'] = m[i]['file2']
		else:
			args[m[i]['library']]['i2'] = tmp['i2']
		if(s_n == 'Demultiplex'):
			if('mid' not in args[m[i]['library']]):
				args[m[i]['library']]['mid'] = {}
			args[m[i]['library']]['mid'][m[i]['SampleID']] = m[i]['mid']
			args[m[i]['library']]['common'] = m[i]['common']
		for k in tmp:
			if(k not in args[m[i]['library']]):
				args[m[i]['library']][k] = tmp[k]
	for library in args:
		args[library]['params'] = p
		args[library]['library'] = library
		_launch_module(args[library],module_name)


def _launch_module(args,module_name):
	module = _create_module(module_name,args)
	if module.execution == 1:
		_exec(module,module_name)
	else:
		log.critical('Skip execution.')

def _exec(module,name):
	if not module.sge:
		for el in module.cmd:
			log.debug(el)
			if log.getLogger().getEffectiveLevel() == 20:
				os.system(el)
	else:
		fw =  open(module.cmd_file, mode='w')
		for el in module.cmd:
			fw.write(el + "\n")
		fw.close()
		qsub_call=''
		if hasattr(module,'iter'):
			if module.iter == 'sample':
				qsub_call = "sbatch -D " + module.wd + " --export=ALL -J " + module.sample + '_' + name + ' --ntasks=' + module.n_cpu + ' ' + module.cmd_file
			elif module.iter == 'library':
				qsub_call = "sbatch -D " + module.wd + " --export=ALL -J " + module.library + '_' + name + ' --ntasks=' + module.n_cpu + ' ' + module.cmd_file
			elif module.iter == 'global':
				qsub_call = "sbatch -D " + module.wd + " --export=ALL -J " + name + ' --ntasks=' + module.n_cpu + ' ' + module.cmd_file
		else:
			if name == 'Blast':
				qsub_call = "sbatch -D " + module.wd + " --export=ALL -J " + module.sample + '_' + name + '_' + module.type + ' ' + module.cmd_file
			else:
				qsub_call = qsub_call = "sbatch -D " + module.wd + " --export=ALL -J " + module.sample + '_' + name + ' ' + ' --ntasks=' + module.n_cpu + ' ' + module.cmd_file
		log.debug(qsub_call)
		if log.getLogger().getEffectiveLevel() == 20:
			os.system(qsub_call)


def _replace_sample_name(step_args,sample_map):
	args = {}
	for k in step_args:
		if(isinstance(step_args[k], str)):
			if(re.search('\(\w+\)',step_args[k])):
				matches = re.findall('\((\w+)\)',step_args[k])
				string = step_args[k]
				for m in matches:
					string = string.replace('('+m+')',sample_map[m])
				args[k] = string
			else:
				args[k] = step_args[k]
		else:
			args[k] = step_args[k]
	return args


def _create_module(name,param):
	if '_' in name:
		name = name.split('_')[0]
	try:
		_class = getattr(importlib.import_module(name),name)
		instance = _class(param)
		return instance
	except Exception as e:
		print("error")
		print(str(e))


def _read_map_file(f):
	reader = csv.reader(f,delimiter="\t")
	data = list(reader)
	headers = data[0]
	headers[0] = headers[0][1:]
	map_obj = []
	for i in range(1,len(data)):
		dict={}
		if len(data[i]) != len(headers):
			print(data[i])
			print(headers)
			sys.exit('line and headers not the same length.')
		for j in range(0,len(headers)):
			dict[headers[j]] = data[i][j]
		map_obj.append(dict)
	return map_obj


def _read_yaml_file(f):
	return yaml.safe_load(f)


def _set_options():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--map',help='The map file.',action='store',type=argparse.FileType('r'),required=True)
	parser.add_argument('-s','--step',help='The step file.',action='store',type=argparse.FileType('r'),required=True)
	parser.add_argument('-n','--name_step',dest='name_step',help='The specified step to launch.',action='store',type=str)
	parser.add_argument('-p','--param',help='The global parameter file.',action='store',type=argparse.FileType('r'))
	parser.add_argument('-v','--verbosity',help='Verbose level', action='store',type=int,choices=[1,2,3,4],default=1)
	args = parser.parse_args()
	return args


if __name__ == "__main__":
	main()
