import csv
import argparse
import logging as log
from collections import defaultdict

def main ():
	options = _set_options()
	_set_log_level(options.verbosity)
	seq_dict = _get_seq(options)
	_print_seq(seq_dict, options.out)


def _get_seq(options):
	seq_dict = defaultdict(list) # dict result
	# for each blastx files
	for i in range(0, len(options.bx_files)):
		blast_file = options.bx_files[i]
		reader = csv.reader(blast_file,delimiter="\t")
		headers = True
		for row in reader:
			# skip headers
			if headers is False:
				contig_name = row[1]
				if row[15] != "":
					seq_dict[contig_name].append(row[15]) # seq
			headers = False
	return(seq_dict)


def _print_seq(sequences, out):
	f = open(out, "w+")
	for seq in sequences:
		f.write(">" + seq + "\n" + '\n'.join(sequences[seq]) + "\n")
	f.close()


def _set_options():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blastx',help='Csv files with blastx results.', action='append',required=True, type=argparse.FileType('r'),dest='bx_files')
	parser.add_argument('-o','--out',help='The output files.',action='store',required=True,type=str,default='./')
	parser.add_argument('-v','--verbosity',help='Verbose level', action='store',type=int,choices=[1,2,3,4],default=1)
	args = parser.parse_args()
	return args

def _set_log_level(verbosity):
	if verbosity == 1:
		log_format = '%(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.INFO,format=log_format)
	elif verbosity == 3:
		log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.DEBUG,format=log_format)

if __name__ == "__main__":
	main()