#!/usr/bin/python3
# -*- coding:utf-8 -*-

current_ver = 'v0.8_220726'

import subprocess
from subprocess import PIPE
import argparse
import os
import tempfile
import time
import sys
import glob # 2020/01/10
import numpy as np
import math
import random
import copy
import csv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Manager
from scipy.spatial import distance


###########################################################################################
###################### Option setting #####################################################
###########################################################################################

def optionswitch():
	parser = argparse.ArgumentParser(description='Make binary sequence from plasmid or geneome sequences')
	parser.add_argument(
		'-g', '--gb',
		type = str,
		dest = 'gb_file_list',
		default = '',
		help = 'GeneBank file list (GB file must be designated with full PATH)'
	) 
	parser.add_argument(
		'-f', '--fasta',
		type = str,
		dest = 'fasta_file_list',
		default = '',
		help = 'Fasta file list (strain<TAB>fasta_file, fasta file must be designated with full PATH)'
	) 
	parser.add_argument(
		'-o', '--out',
		type = str,
		dest = 'outfile',
		default = '',
		help = 'Out put filename'
	) 
	parser.add_argument(
		'--chromosome',
		type = int,
		dest = 'minimum_contig_length',
		default = -1,
		help = 'Make binary for estimated chromosome. Specify the minimum length of contig for analysis. Draft sequence only.'
	) 
	parser.add_argument(
		'--perc_contig',
		type = float,
		dest = 'percent_in_a_contig',
		default = 60,
		help = 'Percent of ORF mapping to chromosome within a contig'
	) 
	parser.add_argument(
		'-e', '--evalue',
		type = float,
		dest = 'evalue',
		default = 0.001,
		help = 'E value (default = 0.001)'
	) 
	parser.add_argument(
		'-c', '--cover',
		type = int,
		dest = 'cover',
		default = 80,
		help = 'Cover ratio against query sequence length (default = 80)'
	) 
	parser.add_argument(
		'-i', '--identity',
		type = int,
		dest = 'sequence_identity',
		default = 80,
		help = 'Threshold of sequence identity (default = 80)'
	) 
	parser.add_argument(
		'--num_process',
		type = int,
		dest = 'the_number_of_process',
		default = 4,
		help = 'The number of process (default = 4)'
	) 
	parser.add_argument(
		'--binary',
		type = str,
		dest = 'prebuild_bin',
		default = '',
	)
	parser.add_argument(
		'--mapping',
		dest = 'mapping',
		action="store_true",
		default = False
	)
	'''
	parser.add_argument(
		'--predict',
		dest = 'predict',
		action="store_true",
		default = False,
		help = 'Predict all ORFs using Prodigal including gb files with annotation data'
	)
	parser.add_argument(
		'--grouping',
		dest = 'grouping',
		action="store_true",
		default = False
	)
	parser.add_argument(
		'--min_share',
		type = float,
		dest = 'shareing_factor',
		default = 50.0,
	)
	parser.add_argument(
		'--j_factor',
		type = float,
		dest = 'joint_factor',
		default = -1.0,
	)
	parser.add_argument(
		'--pickup',
		type = int,
		dest = 'random_pickup',
		default = 0,
	)
	'''
	parser.add_argument(
		'--compatibility',
		dest = 'lower_compatibility',
		action="store_true",
		default = False,
		help = 'Use annotation in GeneBank files or prodigal for files without annotaion.'
	)
	parser.add_argument(
		'-v', '--version',
		dest = 'version_call',
		action="store_true",
		default = False
	)

	args = parser.parse_args()

	gbfiles = args.gb_file_list
	fastafiles = args.fasta_file_list
	out_name = args.outfile
	ev = args.evalue
	min_cover = args.cover
	identity = args.sequence_identity
	#predict_ORF = args.predict
	num_process = args.the_number_of_process
	mapping = args.mapping
	#random_sample = args.random_pickup
	bin_file = args.prebuild_bin
	#g_factor = args.shareing_factor / 100
	#j_factor = args.joint_factor / 100
	#grouping_flag = args.grouping
	contig_length = args.minimum_contig_length
	perc_orf_threshold = args.percent_in_a_contig / 100
	vercall = args.version_call
	lower_compatibility = args.lower_compatibility

	if vercall:
		print(current_ver)
		sys.exit()

	###### option_check #####
	if out_name == '':
		print('-0 or --out is required.')
		sys.exit()

	if gbfiles == '' and fastafiles == '' and bin_file == '':
		print('Designate "gb file list" and/or "fasta file list", or "prebuild binary sequence csv.')
		sys.exit()
	
	return(gbfiles, fastafiles, out_name, ev, min_cover, identity, \
	num_process, mapping, bin_file, \
	contig_length, perc_orf_threshold, \
	vercall, lower_compatibility)
	'''
	return(gbfiles, fastafiles, out_name, ev, min_cover, identity, \
	predict_ORF, num_process, mapping, random_sample, bin_file, \
	g_factor, j_factor, grouping_flag, contig_length, perc_orf_threshold, \
	vercall, lower_compatibility)
	'''
	
###########################################################################################
###################### Read sequences                           ###########################
###########################################################################################

def load_binary(bin_file):
	binary_ndary = np.loadtxt(bin_file, delimiter=',')
	seq_list = []
	info_file = bin_file[0:bin_file.find('_binary.csv')] + '_info.csv'
	with open(info_file, 'r') as infh:
		reader = csv.reader(infh)
		for row in reader:
			seq_list.append(row)
	orfinfo_file = bin_file[0:bin_file.find('_binary.csv')] + '_orfinfo.txt'
	with open(orfinfo_file, 'r') as infh:
		orfinfo = infh.read()
		orf_list = orfinfo.split('\n')
	return(seq_list)

def read_gb(gbfiles):
	acc_num_list =[]
	file = open(gbfiles, "r")
	all_list = file.read().rstrip().split("\n")
	file.close
	print('Reading genbank files')
	strain_temp = []
	remove_tmp = []
	seq_tmp = []
	strain_tmp = []
	while all_list.count("") >0:
		all_list.remove("")
	for each_gb in all_list:
		print(each_gb)
		if os.path.exists(each_gb):
			pass
		else:
			print(each_gb + ' is not found.')
			print(each_gb + ' is removed from the analysis')
			remove_tmp.append(each_gb)
			continue
		record = SeqIO.read(each_gb, "genbank")
		acc_num = record.id.split('.')[0] # Jun. 13, 2018.
		if acc_num in acc_num_list:
			print(acc_num + ' is duplicated.')
			print(each_gb + ' is removed from the analysis')
			remove_tmp.append(each_gb)
			continue
		if len(record.seq) < 20000:
			print(each_gb + ' is smaller than 20000.')
			print(each_gb + ' is removed from the analysis')
			remove_tmp.append(each_gb)
			continue
		if len(record.seq) >= 1000000:
			grouping_flag = False
			
		acc_num_list.append(acc_num)
		seq_tmp.append([record.id, each_gb, 'gb', record.id, 1]) # making seq_list
		strain_tmp.append([record.id, each_gb, 'gb', record.id])
		with open(all_fasta, 'a') as longfasta: # all_fasta containes all sequences as single sequence for one strain
			SeqIO.write(record, longfasta, 'fasta')
			longfasta.flush()
	return(seq_tmp, strain_tmp, remove_tmp)

def read_draft(fastafiles):
	file = open(fastafiles, "r")
	all_draft = file.read().rstrip().split("\n")
	file.close
	print('Reading fasta files')
	strain_temp = []
	while all_draft.count("") >0:
		all_draft.remove("")
	for each_draft in all_draft:
		draft_element = each_draft.split("\t")
		strain_name = draft_element[0]
		draft_path = draft_element[1]
		if os.path.exists(draft_path):
			pass
		else:
			print(draft_path + ' is not found.')
			print(draft_path + ' is removed from the analysis')
			continue
		seq_list.append([strain_name, draft_path, 'fasta', strain_name, 1]) # making seq_list
		strain_name_list.append([strain_name, draft_path, 'fasta', strain_name])
		print(strain_name)
		draft_seqs = SeqIO.parse(draft_path, 'fasta')
		long_seq = ''
		for each_seq in draft_seqs:
			long_seq = long_seq + str(each_seq.seq) + \
			'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn'
		join_draft = SeqRecord(Seq(long_seq), id = strain_name, description = '') # Corrected on 2021/7/30
		with open(all_fasta, 'a') as longfasta:
			SeqIO.write(join_draft, longfasta, 'fasta')
			longfasta.flush()

def read_contig(fastafiles, contig_length):
	file = open(fastafiles, "r")
	all_draft = file.read().rstrip().split("\n")
	file.close
	print('Reading fasta files')
	strain_temp = []
	while all_draft.count("") >0:
		all_draft.remove("")
	for each_draft in all_draft:
		draft_element = each_draft.split("\t")
		strain_name = draft_element[0]
		draft_path = draft_element[1]
		if os.path.exists(draft_path):
			pass
		else:
			print(draft_path + ' is not found.')
			print(draft_path + ' is removed from the analysis')
			continue
		#seq_list.append([seq.id, draft_path, 'fasta', strain_name]) # making seq_list
		print(strain_name)
		strain_name_list.append([strain_name, draft_path, 'fasta', strain_name])
		draft_seqs = SeqIO.parse(draft_path, 'fasta')
		long_seq = ''
		for each_seq in draft_seqs:
			each_seq.id = strain_name + each_seq.id # correct 2021/8/25
			if len(each_seq.seq) >= contig_length:
				seq_list.append([each_seq.id, draft_path, 'fasta', strain_name, 1]) # making seq_list
			else:
				seq_list.append([each_seq.id, draft_path, 'small', strain_name, 0]) # making seq_list
			with open(all_fasta, 'a') as longfasta:
				SeqIO.write(each_seq, longfasta, 'fasta')
				longfasta.flush()



###########################################################################################
###################### genome annotation      2022/7/4         ############################
###########################################################################################
def pickup_orf(strain, seqfile, filetype):
	### for gb file
	if filetype == 'gb':
		record = SeqIO.read(seqfile, "genbank")
		strain = record.id
		temp_fasta = strain + '_temp.fna' 
		with open(temp_fasta, 'w') as gbfasta:
			SeqIO.write(record, gbfasta, 'fasta')
		seqfile = temp_fasta
	
	### dfast annotation
	print(strain)
	dfast_command = ['dfast', '--genome', seqfile, '--out', dfast_dir + strain, '--no_func_anno', '--no_hmm', '--cpu', '1']
	res = subprocess.run(dfast_command)#, stdout=PIPE)

	### pickup ORFs
	dfast_gb = dfast_dir + strain + '/genome.gbk'
	orf_fasta = orf_dir + strain + '_orf.fna' # change dir for orf fasta file 20180123
	with open(orf_fasta, "w") as outfh, open(dfast_gb, 'r') as dfast:
		record = SeqIO.parse(dfast, 'genbank')
		for each_seq in record:
			ft = each_seq.features
			for each_record in ft:
				seq_start = each_record.location.start
				seq_end = each_record.location.end
				length = seq_end - seq_start
				if each_record.type == "CDS" and length <= 10000:
					sub_gb_record = each_seq[seq_start:seq_end]
					sub_gb_record.description = sub_gb_record.id
					sub_gb_record.id = strain + '_' + each_record.qualifiers.get('locus_tag',"n")[0] +\
					"-" + str(seq_start) + "-" + str(seq_end) # ORF naming manner : seq.id + ORF_number
					SeqIO.write(sub_gb_record, outfh, "fasta")
	if filetype == 'gb':
		os.remove(temp_fasta)
	return(orf_fasta)

def gb_orf(gb):
	record = SeqIO.read(gb, "genbank")
	print(record.id)
	ft = record.features
	orf_fasta = orf_dir + record.id + '_orf.fna' # change dir for orf fasta file 20180123
	if len(ft) <= 2: # or predict_ORF:
		no_annotation = open(error_file, 'a') 
		no_annotation.write(gb)
		no_annotation.write('\n')
		no_annotation.close
		record = SeqIO.read(gb, "genbank")
		gb_fasta = orf_dir + record.id + ".fasta" # 20160125
		with open(gb_fasta, 'w') as gbfasta:
			SeqIO.write(record, gbfasta, 'fasta')
		prodigal_command = ['prodigal', '-i', gb_fasta, '-d', orf_fasta, '-q']
		res = subprocess.run(prodigal_command, stdout=PIPE)
		return(orf_fasta)
	with open(orf_fasta, "w") as outfh:
		for each_record in ft:
			seq_start = each_record.location.start
			seq_end = each_record.location.end
			length = seq_end - seq_start
			if each_record.type == "CDS" and length <= 10000:
				sub_gb_record = record[seq_start:seq_end]
				sub_gb_record.description = sub_gb_record.id
				sub_gb_record.id = record.id + '_' + each_record.qualifiers.get('locus_tag',"n")[0] +\
				"-" + str(seq_start) + "-" + str(seq_end)
				SeqIO.write(sub_gb_record, outfh, "fasta")
	return(orf_fasta)

def draft_orf(strain, draft):
	print(draft)
	orf_fasta = orf_dir + strain + '_orf.fna' # change dir for orf fasta file 20180123
	prodigal_command = ['prodigal', '-i', draft, '-d', strain + '_temp.fna', '-q']
	res = subprocess.run(prodigal_command)#, stdout=PIPE)
	with open(orf_fasta, 'w') as outfasta, open(strain + '_temp.fna', 'r') as infasta:
		seq = SeqIO.parse(infasta, "fasta")
		for each_seq in seq:
			each_seq.id = strain + "_" + each_seq.id
			SeqIO.write(each_seq, outfasta, "fasta")
	os.remove(strain + '_temp.fna')
	return(orf_fasta)



###########################################################################################
###################### ORF check                                            ###############
###################### These modules work only in this program              ###############
###########################################################################################

def orf_check(orf_fasta, nonredundant): # remove redundant ORFs, 20180123
	blast_temp = tempfile.NamedTemporaryFile()
	db_name = blast_temp.name
	blast_temp.close
	mkbl_command = ['makeblastdb', '-dbtype', 'nucl', '-hash_index', \
	'-in', orf_fasta, '-out', db_name]
	res = subprocess.run(mkbl_command, stdout=PIPE)
	blastn_cline = NcbiblastnCommandline(query=orf_fasta, db=db_name, evalue=ev,\
	outfmt=6, perc_identity=identity, qcov_hsp_perc=min_cover)
	stdout, stderr = blastn_cline()
	blast_results = stdout.rstrip().split('\n')
	stdout=''
	dbfiles = glob.glob(db_name + '.n*') #2020/01/10
	for dbfile in dbfiles:
		os.remove(dbfile)
	
	orf_tempdic = {}
	i = 0
	with open(orf_fasta, 'r') as orf:
		sequence = SeqIO.parse(orf, 'fasta')
		for each in sequence:
			orf_tempdic[each.id] = [0, each.id, each.description, each.seq, i]
			i += 1
	for each_r in blast_results:
		br = each_r.split('\t')
		if orf_tempdic[br[0]][0] >= 0:
			orf_tempdic[br[0]][0] = 1
			if br[0] != br[1]:
				orf_tempdic[br[1]][0] = -1
	with open(nonredundant, 'w') as newfasta:
		for each_seq, values in sorted(orf_tempdic.items(), key=lambda x: x[1][4]):
			if values[0] == 1:
				seq_items = SeqRecord(values[3], id = values[1], description = values[2])
				SeqIO.write(seq_items, newfasta, 'fasta')

def pickup_unique(strain_id):
	orf_tempdic = {}

	if lower_compatibility:
		strain = strain_name_list[strain_id][0]
		if strain_name_list[strain_id][2] == 'gb': 
			orf_fasta = gb_orf(strain_name_list[strain_id][1]) # recieve orf_fasta name from gb_orf() 20210817
		elif strain_name_list[strain_id][2] == 'fasta':
			orf_fasta = draft_orf(strain, strain_name_list[strain_id][1]) # recieve orf_fasta name from draft_orf() 20210817
	else:
		orf_fasta = pickup_orf(strain_name_list[strain_id][0], strain_name_list[strain_id][1], strain_name_list[strain_id][2])
	
	nonredundant_fasta = orf_dir + strain_name_list[strain_id][0] + '_nonredundant.fna' # change dir of orf fasta file 20180123
	orf_check(orf_fasta, nonredundant_fasta)
	nonredundant_num = strain_id * 100000 # orf id
	
	if strain_id == 0: #20180123
		with open(nonredundant_fasta, 'r') as orf:
			sequence = SeqIO.parse(orf, 'fasta')
			for each in sequence:
				orf_dic[nonredundant_num] = [each.id, each.description, each.seq]
				nonredundant_num += 1
		return

	with open(nonredundant_fasta, 'r') as orf:
		sequence = SeqIO.parse(orf, 'fasta')
		for each in sequence:
			orf_tempdic[each.id] = [0, each.id, each.description, each.seq]
	
	#seq_list [seq.id, path_to_file, 'gb/fasta', strain_name]
	#strain_name_list [strain_name, path_to_file, 'gb/fasta', strain_name]
	seqid_temp = tempfile.NamedTemporaryFile()
	with open(seqid_temp.name, 'w') as use_seq:
		seq_temp = []
		for i in range(strain_id):
			seq_temp.append(strain_name_list[i][0]) # append strain_name
		for each in seq_list: # make list adapting to blast search 180123
			if each[3] in seq_temp:
				use_seq.write(each[0] + "\n")  # append seq.id to seqidlist for blastn
		
	blastn_cline = NcbiblastnCommandline(query=nonredundant_fasta, db=out_name, evalue=ev,\
	outfmt=6, seqidlist=seqid_temp.name, perc_identity=identity, qcov_hsp_perc=min_cover,\
	num_alignments=2) # blast db "out_name" is constructed from _all.fasta
	stdout, stderr = blastn_cline()
	if stdout == '':
		for each_seq, values in orf_tempdic.items(): # orf_tempdic[each.id] = [0, each.id, each.description, each.seq]
			orf_dic[nonredundant_num] = [values[1], values[2], values[3]]
			nonredundant_num += 1
		return

	blast_results = stdout.rstrip().split('\n')
	stdout=''
	seqid_temp.close()

	for each_r in blast_results:
		br = each_r.split('\t')
		orf_tempdic[br[0]][0] += 1
	for each_seq, values in orf_tempdic.items():
		if values[0] == 0:
			orf_dic[nonredundant_num] = [values[1], values[2], values[3]]
			nonredundant_num += 1

def run_pickup_unique(all_fasta, out_name, orf_dic):
	mkbl_command = ['makeblastdb', '-dbtype', 'nucl', '-parse_seqids', \
	'-in', all_fasta, '-out', out_name]
	res = subprocess.run(mkbl_command, stdout=PIPE)
	orf_list = []
	orf_path = out_name + '_orf'
	if os.path.isdir(orf_path): #20180123
		print('Create fasta files in "orf" directory')
	else:
		print('Create "orf" directory')
		os.mkdir(orf_path)
		print('Create fasta files in "orf" directory')
	#orf_dir = current_dir + '/' + orf_path + '/'
	pool = Pool(num_process)	
	pool.map(pickup_unique, range(num_strain))
	pool.close()
	num_orf = len(orf_dic)
	return(orf_list, num_orf)

def write_elemetORF(elementORF_fasta, orf_dic, query_temp):
	print('Writing elemental ORFs') # 20180125
	#query_temp = tempfile.NamedTemporaryFile() # Warning!!! bad use of query_temp
	with open(elementORF_fasta, 'w') as outfh, open(query_temp.name, 'w') as queryfh: # Warning!!! bad use of query_temp
		i = 0
		for key, values in sorted(orf_dic.items(), key=lambda x: x[0]):
			seq_items = SeqRecord(values[2], id = values[0], description = values[1])
			SeqIO.write(seq_items, outfh, 'fasta')
			query_items = SeqRecord(values[2], id = str(i)) # Sequences are renamed with serial number
			SeqIO.write(query_items, queryfh, 'fasta')
			orf_list.append(values[0])
			i += 1

###########################################################################################
###################### binarize ###########################################################
###########################################################################################

def make_binary(seq_id):
	#seq_list [seq.id, path_to_file, 'gb/fasta', strain_name, adapt]
	#strain_name_list [strain_name, path_to_file, 'gb/fasta', strain_name]
	print(seq_list[seq_id][0])
	seqid_temp = tempfile.NamedTemporaryFile()
	with open(seqid_temp.name, 'w') as use_seq:
		use_seq.write(seq_list[seq_id][0])
	blastn_cline = NcbiblastnCommandline(query=query_temp.name, db=out_name, evalue=ev,\
	outfmt=6, seqidlist=seqid_temp.name, perc_identity=identity, qcov_hsp_perc=min_cover,\
	num_alignments=1)
	try:
		stdout, stderr = blastn_cline()
	except Exception as e:
		pass
	else:
		if stdout == '':
			return()
		blast_results = stdout.rstrip().split('\n')
		stdout=''
		seqid_temp.close()
		for each_r in blast_results:
			br = each_r.split('\t')
			positive_list.append([seq_id, int(br[0])])	


###########################################################################################
###################### Write nexus file ###################################################
###########################################################################################
'''
def group_nex(fnum, glist, group_factor, join_factor, random_sample_nex):
	if fnum == 999:
		g_nex_filename = random_sample_nex + '.nex'
		g_fasta_filename = random_sample_nex + '.fasta'
	else:
		g_nex_filename = out_name + '_group' + str(fnum).zfill(3) + '.nex' 
		g_fasta_filename = out_name + '_group' + str(fnum).zfill(3) + '.nex' 
	temp_binary = binary_ndary[glist,:].astype(np.int)  ###### binary of selected strains
	orf_excluded_num = np.sum(temp_binary, axis = 0).astype(np.int)
	exclude_orf = []
	for i in range(orf_excluded_num.shape[0]):
		if orf_excluded_num[i] == 0:
			exclude_orf.append(i)
	temp_binary = np.delete(temp_binary, exclude_orf, 1)
	print('Writing binary sequence : ' + g_nex_filename)
	nexus(g_nex_filename, g_fasta_filename, temp_binary, strain_list, glist, 1, group_factor, join_factor)
'''
def nexus(out_file, binfasta_name, binaryseqs, seq_list, write_seq_list, group_flag, group_factor, join_factor):
	
	strain = len(write_seq_list) #binaryseqs.shape[0]
	orf = binaryseqs.shape[1]
	div_lines = int(math.ceil(orf / 60.0))

	contents = open(out_file, 'w')
	fasta = open(binfasta_name, 'w')
	contents.write('#nexus\n')
	contents.write("[!" + out_file + "]\n")
	contents.write("begin taxa;\n")
	contents.write("dimensions ntax=" + str(strain) + ";\n")
	contents.write("end;\n")
	contents.write("\n")
	contents.write("begin characters;\n")
	contents.write("dimensions nchar=" + str(orf) + ";\n")
	contents.write("format datatype=standard labels;\n")
	contents.write("matrix\n")
	for genome_id in range(strain):
		contents.write(seq_list[write_seq_list[genome_id]][0] + "\n")
		fasta.write('>' + seq_list[write_seq_list[genome_id]][0] + "\n")
		for i in range(div_lines): #20180129 ##textwrap is very slow in cases of large sequences
			sequence = ""
			seq_start = i * 60
			seq_end = seq_start + 60
			if seq_end >= orf:
				seq_end = orf
			for each_binary in range(seq_start, seq_end):
				sequence = sequence + str(int(binaryseqs[genome_id, each_binary]))
			contents.write(sequence + "\n")
			fasta.write(sequence + "\n")
	contents.write(";\n")
	contents.write("end;\n")
	contents.write("\n")
	contents.write("[!\n") #begin comments
	contents.write("------- " + current_ver + " ----------------------------------------\n")
	contents.write("------- blastn settings ------------------------------------\n")
	contents.write("evalue = " + str(ev) + "\n")
	contents.write("perc_identity = " + str(identity) + " %\n")
	contents.write("qcov_hsp_perc = " + str(min_cover) + " %\n")
	contents.write("\n")
	if group_flag:
		contents.write("------- Grouping option ------------------------------------\n")
		contents.write("minimum share percentage (--min_share) : " + str(group_factor * 100) + " %\n")
		contents.write("minimum percentage to connect subgroups (--j_factor) : " + str(join_factor * 100) + " %\n")
		contents.write("\n")
	contents.write("------- Sequence files -------------------------------------\n")
	for eachfile in write_seq_list:
		contents.write(seq_list[eachfile][0] + "\t" + seq_list[eachfile][1] + "\n")
	contents.write("------------------------------------------------------------\n")
	contents.write("]\n") # end comments

	
	contents.flush()
	contents.close
	fasta.flush()
	fasta.close()


###########################################################################################
###################### Statistics #########################################################
###########################################################################################

def calc_homology(seq1):
	second_seq = seq1 + 1
	dice_temp_list = [1 - distance.dice(binary_ndary[seq1],binary_ndary[i]) for i in range(num_sequence)]
	dice_dic[seq1] = dice_temp_list
	share_temp_list = [np.sum(np.logical_and(binary_ndary[seq1], binary_ndary[i])) for i in range(num_sequence)]
	share_dic[seq1] = share_temp_list

def part_share_indices(index_list):
	share_ndary = share_indices_ndary[index_list,:]
	share_ndary = share_ndary[:,index_list]
	return(share_ndary)

'''
###########################################################################################
###################### Grouping ###########################################################
###########################################################################################

def grouping(sf, jf):
	remain_seq_list = list(range(num_sequence))  #modify 2020/5/1
	remain_seq_2nd = list(range(num_sequence))   #modify 2020/5/1
	group_2nd = []
	for i in range(len(remain_seq_list)): #pick up closely related sequences
		temp_group = []
		for j in range(i + 1, len(remain_seq_list)):
			smaller_seq = min(share_indices_ndary[remain_seq_list[i]][remain_seq_list[i]], \
			share_indices_ndary[remain_seq_list[j]][remain_seq_list[j]])
			if share_indices_ndary[remain_seq_list[i]][remain_seq_list[j]] >= smaller_seq * sf:
				temp_group.append(remain_seq_list[j])
				if remain_seq_list[j] in remain_seq_2nd:
					remain_seq_2nd.remove(remain_seq_list[j])
		if len(temp_group) >= 1:
			temp_group.append(remain_seq_list[i])
			group_2nd.append(temp_group)
			if remain_seq_list[i] in remain_seq_2nd:
				remain_seq_2nd.remove(remain_seq_list[i])
	pre_len_group_2nd = len (group_2nd)
	post_len_group_2nd = 0
	while pre_len_group_2nd > post_len_group_2nd:
		pre_len_group_2nd = len (group_2nd)
		group_2nd = combine_group(group_2nd, jf)
		post_len_group_2nd = len(group_2nd)
	
	return(group_2nd, remain_seq_2nd)

	
def combine_group(g_list, parameter1):
	# combine similar groups
	len_g_list = len(g_list) # 
	for i in range(len_g_list - 1): #Compare 1st and 2nd groups one by one
		if len(g_list[i]) == 0: #Skip empty group
			pass
		else:
			for j in range(i + 1, len_g_list):
				k = 0
				len_smaller_group = min([len(g_list[i]), len(g_list[j])])		
				for each in g_list[j]: #Count the number of shared seqs in 2nd group
					if each in g_list[i]:
						k += 1
				if k >=  len_smaller_group * parameter1: # If more than X% of indices in 2nd group share 1st   
					g_list[i].extend(g_list[j]) #Combine 1st and 2nd group
					g_list[i] = list(set(g_list[i])) #Remove overlapped indices
					g_list[j] = [] # Remove 2nd group
	# Remove empty group
	del_list = []
	for i in range(len(g_list)):
		if len(g_list[i]) == 0:
			del_list.append(i)
	del_list.reverse()
	for i in del_list:
		del g_list[i]
	return(g_list)
'''

###########################################################################################
###################### Write mapping file #################################################
###########################################################################################

def write_map(out_name, binary_ndary, seq_list, orf_list, mask):
	num_sequence = binary_ndary.shape[0]
	num_orf = binary_ndary.shape[1]
	mask = mask.tolist()
	contents = open(out_name, 'w')
	contents.write("\t")
	for i in range(num_sequence):
		contents.write(seq_list[i][0] + "\t")
	contents.write("sum\t")
	contents.write("chromosome\n")
	for i in range(num_orf):
		orf_sum = 0
		contents.write(orf_list[i] + "\t")
		for j in range(num_sequence):
			contents.write(str(binary_ndary[j, i]) + "\t")
			orf_sum = orf_sum + binary_ndary[j, i]
		contents.write(str(orf_sum)+ '\t')
		contents.write(str(mask[i])+ '\n')
	contents.flush()
	contents.close

'''
def out_mapping(group_list, out_fh):
	print('Making mapping file')
	temp_ndary = binary_ndary[group_list]
	temp_num_orf = temp_ndary.shape[1]
	temp_num_strain = temp_ndary.shape[0]
	contents = open(out_fh + '_mapping.txt', 'w')
	contents.write("\t")
	for genome_id in group_list:
		contents.write(seq_list[genome_id][0] + "\t")
	contents.write("sum\n")
	for each_binary in range(temp_num_orf):
		orf_sum = 0
		contents.write(orf_list[each_binary] + "\t")
		for genome_id in range(temp_num_strain):
			contents.write(str(temp_ndary[genome_id, each_binary]) + "\t")
			orf_sum = orf_sum + temp_ndary[genome_id, each_binary]
		contents.write(str(orf_sum)+ '\n')
	contents.flush()
	contents.close
	print('Done! Mapping file : ' + out_fh + '_mapping.txt')
'''
def write_list(outfile, in_list, header):
	with open(outfile, 'w') as outf:
		writer = csv.writer(outf, delimiter = '\t')
		if header != '':
			writer.writerow(header)
		writer.writerows(in_list)

def write_line(outfile, in_list):
	with open(outfile, 'w') as outfh:
		for each in in_list:
			outfh.write(each + '\n')


###########################################################################################
###################### Chromosome typing 2021/12/3 ########################################
###########################################################################################

def check_chromosome(num_sequence, seq_list):
	gb_list = []
	contig_list = []
	for i in range(num_sequence):
		if seq_list[i][2] in ['gb', 'chromosome']:
			gb_list.append(i)
		else:
			contig_list.append(i)
	return(gb_list, contig_list)
	# gb_list: sequence number of chromosome
	# contig_list: sequence number of unknown contigs

def make_chromosome_mask(gb_list, binary_ndary, modified_seq_list, num_orf):
	chr_mask_ndary = np.zeros(num_orf, dtype = int)
	for i in gb_list:
		chr_mask_ndary = np.logical_or(chr_mask_ndary, binary_ndary[i:i+1])
		modified_seq_list[i][2] = 'chromosome'
	return(modified_seq_list, chr_mask_ndary)
	# modified_seq_list: seq_list with chromosome info [strain, path, 'chromosome', strain/AccNo, 1]
	# chr_mask_ndary: 1D ndarray put '1' for chromosomal ORFs

def judge_chromosome(contig_list, binary_ndary,chr_mask_ndary, modified_seq_list, threshold):
	for i in contig_list:
		orf_sum = np.sum(binary_ndary[i:i+1])
		masked_sum = np.sum(binary_ndary[i:i+1][chr_mask_ndary])
		modified_seq_list[i].append(int(orf_sum))
		modified_seq_list[i].append(int(masked_sum))
		if orf_sum == 0:
			modified_seq_list[i][2] = 'none'
		elif masked_sum / orf_sum >= threshold and modified_seq_list[i][4] == 1:
			modified_seq_list[i][2] = 'chromosome'
		elif modified_seq_list[i][4] == 1:
			modified_seq_list[i][2] = 'plasmid'
		else:
			modified_seq_list[i][2] = 'none'
	return(modified_seq_list)
	# modified_seq_list: [strain, path, "updated info", strain/AccNo, 1, orf_sum, masked_sum]

def make_new_seq_list(modified_seq_list):
	new_seq_list = []
	for each in modified_seq_list:
		if each[3] not in new_seq_list:
			new_seq_list.append(each[3])
	for i in range(len(new_seq_list)):
		new_seq_list[i] = [new_seq_list[i], '', [], []]

	j = 0
	for i in range(len(modified_seq_list)):
		if modified_seq_list[i][3] != new_seq_list[j][0]:
			j += 1
		if modified_seq_list[i][2] == 'chromosome':
			new_seq_list[j][2].append(i)
		elif modified_seq_list[i][2] == 'plasmid':
			new_seq_list[j][3].append(i)
		new_seq_list[j][1] = modified_seq_list[i][1]
	return(new_seq_list)
	# new_seq_list: [strain/AccNo, path, [i1, i2, ... for chromosome], [i1, i2, ... for plasmid]]

def chr_binary(binary_ndary, new_seq_list):
	num_orf = binary_ndary.shape[1]
	binary_chr_ndary = np.zeros((len(new_seq_list), num_orf), dtype = int) # 2D ndarray
	binary_pls_ndary = np.zeros((len(new_seq_list), num_orf), dtype = int) # 2D ndarray
	binary_whole_ndary = np.zeros((len(new_seq_list), num_orf), dtype = int) # 2D ndarray
	chromosome_mask = np.zeros(num_orf, dtype = int) # 1D ndarray
	plasmid_mask = np.zeros(num_orf, dtype = int) # 1D ndarray

	for i in range(len(new_seq_list)): # i: strain No.
		for j in new_seq_list[i][2]: # j: each contig No. "chromossome"
			#binary_chr_ndary[i:i+1] = np.logical_or(binary_chr_ndary[i:i+1], binary_ndary[j:j+1]) # pile up binary of each contig
			chromosome_mask[0:num_orf] = np.logical_or(chromosome_mask[0:num_orf], binary_ndary[j:j+1]) # pile up for chromosome_mask
		for j in new_seq_list[i][3]: # j: each contig No. "plasmid"
			#binary_pls_ndary[i:i+1] = np.logical_or(binary_pls_ndary[i:i+1], binary_ndary[j:j+1]) # pile up binary of each contig
			plasmid_mask[0:num_orf] = np.logical_or(plasmid_mask[0:num_orf], binary_ndary[j:j+1])  # pile up for plasmid_mask
		whole_list = new_seq_list[i][2] + new_seq_list[i][3]
		for j in whole_list:
			binary_whole_ndary[i:i+1] = np.logical_or(binary_whole_ndary[i:i+1], binary_ndary[j:j+1]) # pile up binary of all contigs
		binary_pls_ndary = np.logical_and(binary_whole_ndary, plasmid_mask).astype(np.int)
		non_plasmid_mask = np.logical_not(plasmid_mask) # make mask to exclude ORFs found in plasmids
		#chromosome_mask = np.logical_and(chromosome_mask, non_plasmid_mask) # modify chromosome_mask to exclude ORFs found in plasmids
		binary_chr_ndary = np.logical_and(binary_whole_ndary, chromosome_mask).astype(np.int) # pick up ORFs found in chromosome
		binary_chr_ndary = np.logical_and(binary_chr_ndary, non_plasmid_mask).astype(np.int) # exclude ORFs found in plasmids
		# "Chromosome" binary is defied by chromosome mask to avoid assembly errors and recombinations 21/12/3
		chromosome_mask = np.logical_and(chromosome_mask, non_plasmid_mask) # 2022/7/26
	return(binary_chr_ndary, binary_pls_ndary, binary_whole_ndary, chromosome_mask)

def chr_typing(seq_list, num_sequence, binary_ndary, perc_orf_threshold, out_name, orf_list):
	modified_seq_list = copy.copy(seq_list) # [seq.id, draft_path, 'gb/fasta/small', strain_name]
	gb_list, contig_list = check_chromosome(num_sequence, seq_list)
	modified_seq_list, chr_mask_ndary = make_chromosome_mask(gb_list, binary_ndary, modified_seq_list, num_orf)
	modified_seq_list = judge_chromosome(contig_list, binary_ndary, chr_mask_ndary, modified_seq_list, perc_orf_threshold)
	
	################ Make new sequence list ###############
	new_seq_list = make_new_seq_list(modified_seq_list)

	################ Make binary for "chromosome" and "plasmid" each ###############
	binary_chr_ndary, binary_pls_ndary, binary_whole_ndary, mask_ndary = chr_binary(binary_ndary, new_seq_list)
	num_sequence = len(new_seq_list)
	header_tup = ('contig', 'file', 'replicon', 'strain', 'use', 'Total ORF', 'Chromosome ORF')
	write_list(out_name + '_contig_info.txt', modified_seq_list, header_tup)
	nexus(out_name + '_chromosome.nex', out_name + '_chromosome_bin.fasta', \
	binary_chr_ndary, new_seq_list, list(range(num_sequence)), 0, 0, 0)
	nexus(out_name + '_plasmid.nex', out_name + '_plasmid_bin.fasta', \
	binary_pls_ndary, new_seq_list, list(range(num_sequence)), 0, 0, 0)
	if mapping:
		print('Making mapping file')
		#write_map(out_name + '_all_mapping.txt', binary_ndary, seq_list, orf_list, mask_ndary)
		write_map(out_name + '_chromosome_mapping.txt', binary_chr_ndary, new_seq_list, orf_list, mask_ndary)
		write_map(out_name + '_plasmid_mapping.txt', binary_pls_ndary, new_seq_list, orf_list, mask_ndary)
	return(binary_whole_ndary, mask_ndary, num_sequence, new_seq_list)


#####################################################################################################
#####   blast output information 
#####   br[0] query acc.ver
#####   br[1] subject acc.ver
#####   br[2] % identity
#####   br[3] alignment length
#####   br[4] mismatches
#####   br[5] gap opens
#####   br[6] q. start
#####   br[7] q. end
#####   br[8] s. start
#####   br[9] s. end
#####   br[10] evalue
#####   br[11] bit score
#####################################################################################################






if __name__ == '__main__':
	
	############# initialize ###########################################
	gbfiles, fastafiles, out_name, ev, min_cover, identity, \
	num_process, mapping, bin_file, \
	contig_length, \
	perc_orf_threshold, vercall, lower_compatibility = optionswitch()
	'''
	gbfiles, fastafiles, out_name, ev, min_cover, identity, \
	predict_ORF, num_process, mapping, random_sample, bin_file, \
	g_factor, j_factor, grouping_flag, contig_length, \
	perc_orf_threshold, vercall, lower_compatibility = optionswitch()
	'''
	start = time.time()
	elementORF_fasta = out_name + "_elementORF.fasta"
	all_fasta = out_name + '_all.fasta'
	error_file = out_name + "_no_annotation.txt"
	current_dir = os.getcwd() # set current directory 20180123
	
	seq_dict = {}
	seq_list = [] # [seq.id, draft_path, 'gb/fasta/small', strain_name, 1/0] (0 for small contig)
	strain_name_list = []
	manager = Manager()
	orf_dic = manager.dict() # unique orf contaner 20180124
	
	############# load files ###########################################
	if bin_file != '':
		seq_list = load_binary(bin_file)

	if gbfiles != '' and bin_file == '':
		seq_list, strain_name_list, remove_seq_list = read_gb(gbfiles)
		
	if fastafiles != "" and contig_length < 0 and bin_file == '':
		read_draft(fastafiles)
	
	if contig_length >= 0:
		read_contig(fastafiles, contig_length)

	num_sequence = len(seq_list)
	num_strain = len(strain_name_list)

	############# pick up unique ORFs ##################################
	if bin_file == '':
		print('file summary writing')
		orf_dir = current_dir + '/' + out_name + '_orf' + '/'
		dfast_dir = current_dir + '/' + out_name + '_dfast' + '/'
		orf_list, num_orf = run_pickup_unique(all_fasta, out_name, orf_dic)
		query_temp = tempfile.NamedTemporaryFile() # Warning!!! bad use of query_temp
		write_elemetORF(elementORF_fasta, orf_dic, query_temp)
	

		############# making binary sequences #################################
		print('Making binary sequences')
		positive_list = manager.list()
		pool = Pool(num_process)	
		pool.map(make_binary, range(num_sequence))
		pool.close()
		query_temp.close()                        # Warning!!! bad use of query_temp
		print('Constructiong binary array')
		binary_ndary = np.zeros((num_sequence, num_orf), dtype = int) ##### create binary sequence array
		for values in positive_list:
			binary_ndary[values[0]][values[1]] = 1

		dbfiles = glob.glob(out_name + '*.n*') #2020/01/10
		for dbfile in dbfiles:
			os.remove(dbfile)

	orf_dic.clear

	############# Binary for chromosome ####################################
	if contig_length >= 0:
		binary_whole_ndary, mask_ndary, num_sequence, seq_list = \
		chr_typing(seq_list, num_sequence, binary_ndary, perc_orf_threshold, out_name, orf_list)
		binary_ndary = binary_whole_ndary

	
	############# writing mapping file ####################################
	if mapping and bin_file == '':
		print('Making mapping file')
		mask_ndary = np.ones(num_orf, dtype = int)
		write_map(out_name + '_all_mapping.txt', binary_ndary, seq_list, orf_list, mask_ndary)
		print('Done!')

	############# Calculating homology #############################
	print('Calculating homologies')
	share_dic = manager.dict()
	dice_dic = manager.dict()
	pool = Pool(num_process)	
	pool.map(calc_homology, range(num_sequence))
	pool.close()
		
	share_indices_ndary = np.empty((0, num_sequence), int)
	index_no = 0
	for i in range(num_sequence):
		share_indices_ndary = np.append(share_indices_ndary, [share_dic[i]], axis = 0)

	'''
	############# Making group #####################################
	if grouping_flag:
		if j_factor < 0:
			j_factor = 0.09
			i = 1
			threshold = 0
			while i > 0:
				i = 0
				j_factor += 0.01
				if j_factor >= 0.5:
					threshold = larger_num
					print('threshold is changed to ' + str(threshold))
					j_factor = 0.1
				print('--j_factor = ' + str(j_factor * 100))
				group_list, remain_seq = grouping(g_factor, j_factor)
				larger_num = 0
				for each_g in group_list:
					if len(each_g) >= num_sequence * 0.25: # evaluate only groups larger then num_sequence * 0.25
						temp_share_index = part_share_indices(each_g) # make array including a group
						share_median_ndary = np.median(temp_share_index, axis = 0) # make array of median value of share_index for each strain
						num_lowsim = np.sum(share_median_ndary < 10)
						if num_lowsim > threshold:
							larger_num = max([num_lowsim, larger_num])
							print(larger_num, num_lowsim)
							i += 1
		else:
			group_list, remain_seq = grouping(g_factor, j_factor)

		i = 0
		len_group_list = []
		for each_g in group_list:
			len_group_list.append([len(each_g), i])
			i += 1
		len_group_list.sort()
		len_group_list.reverse()
		filenum = 0
		for each_group in len_group_list:
			out_head = out_name + '_group' + str(filenum).zfill(3)
			list_file = out_head + '.list'
			print(list_file)
			with open(list_file, 'w') as outfh:
				for seq in group_list[each_group[1]]:
					outfh.write(seq_list[seq][1] + '\n')
			if len(group_list[each_group[1]]) >= 3: #and len(group_list) >= 2:
				group_nex(filenum, group_list[each_group[1]], g_factor, j_factor)
				if mapping:
					out_mapping(group_list[each_group[1]], out_head)
				temp_share_index = part_share_indices(group_list[each_group[1]])
				out_csv = out_head + '_share_index.csv'
				print(out_csv)
				np.savetxt(out_csv, temp_share_index, delimiter=',')
			filenum += 1

		minor_out = out_name + '_minor.txt'
		print('Minor sequences are listed in ' + minor_out)
		with open(minor_out, 'w') as outfh:
			for each in remain_seq:
				outfh.write(seq_list[each][0] + '\t' + seq_list[each][1] + '\n')
		
		############# random select ############################################
		if random_sample > 0:
			i = 0
			for seq in group_list:
				if len(group_list[i]) > random_sample:
					random_sample_nex = out_name  + '_group' + str(i).zfill(3) + \
					'_random' + str(random_sample)
					randomselect = random.sample(group_list[i], random_sample)
					group_nex(999, randomselect, g_factor, j_factor, random_sample_nex)#correct 181115
					with open(random_sample_nex + '.nex', 'a') as outfh:
						outfh.write("[!\n")
						outfh.write("------- random sampling ------------------------------------\n")
						outfh.write(str(random_sample) + " sequences are selected randomly from group" + \
						str(i).zfill(3) + "\n")
						outfh.write("]\n")
				i += 1
	'''
	############# Writing whole binary sequence info #######################
	if bin_file == '':
		############# Output Nexus file ###########################
		seq_num_list = range(num_sequence)
		print('Writing whole binary sequence')
		out_nexus = out_name + '_all.nex'
		binfasta_name = out_name + '_all_bin.fasta'
		nexus(out_nexus, binfasta_name, binary_ndary, seq_list, seq_num_list, 0, 0, 0)
		print('Done! Binary sequence file : ' + out_nexus)

		if bin_file == '':
			np.savetxt(out_name + '_binary.csv', binary_ndary, delimiter=',')
			write_list(out_name + '_info.csv', seq_list, '')
			write_line(out_name + '_orfinfo.txt', orf_list)

		############# Output Dice indices #########################
		size = binary_ndary.shape
		seq_length = size[1]
		dice_line = ''
		for seq_num in range(len(seq_list)):
			dice_line = dice_line  + '\t'+ seq_list[seq_num][0]
		dice_list = [dice_line]
		d_indices = np.array([])
		for i in range(num_sequence):
			mapped_dice = map(str, dice_dic[i])
			dice_line = '\t'.join(mapped_dice)
			dice_line = seq_list[i][0] + '\t' + dice_line
			dice_list.append(dice_line)
			d_indices = np.append(d_indices, np.array(dice_dic[i])) 
		write_line(out_name + "_dice.txt", dice_list)
		print('Done! Dice indeces file : ' + out_name + '_dice.txt')

		############# Output Share indices #########################
		share_line = ''
		for seq_num in range(len(seq_list)):
			share_line = share_line  + '\t'+ seq_list[seq_num][0]
		share_list = [share_line]
		for i in range(num_sequence):
			mapped_share = map(str, share_dic[i])
			share_line = '\t'.join(mapped_share)
			share_line = seq_list[i][0] + '\t' + share_line
			share_list.append(share_line)
		write_line(out_name + "_share.txt", share_list)
		print('Done! Share indeces file : ' + out_name + '_share.txt')

	############# Closing ##################################################
	elapsed_time = time.time() - start
	print("elapsed_time:{0}".format(elapsed_time) + "[sec]")
