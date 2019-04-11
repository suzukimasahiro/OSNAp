#!/usr/bin/python
# -*- coding:utf-8 -*-

import subprocess
import argparse
import os
import tempfile
import time
import sys
import numpy as np
import math
import random
import copy
import csv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Manager
from scipy.spatial import distance

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
	required = True, 
	help = 'Out put filename'
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
	'--predict',
	dest = 'predict',
	action="store_true",
	default = False,
	help = 'Predict all ORFs using Prodigal including gb files with annotation data'
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

args = parser.parse_args()

gbfiles = args.gb_file_list
fastafiles = args.fasta_file_list
out_name = args.outfile
ev = args.evalue
min_cover = args.cover
identity = args.sequence_identity
predict_ORF = args.predict
num_process = args.the_number_of_process
mapping = args.mapping
random_sample = args.random_pickup
bin_file = args.prebuild_bin
g_factor = args.shareing_factor / 100
j_factor = args.joint_factor / 100
grouping_flag = args.grouping

current_ver = 'v0.6_190225'

def gb_orf(gb):
	record = SeqIO.read(gb, "genbank")
	print(record.id)
	ft = record.features
	orf_fasta = orf_dir + record.id + '_orf.fna' # change dir for orf fasta file 20180123
	if len(ft) <= 2 or predict_ORF:
		no_annotation = open(error_file, 'a') 
		no_annotation.write(gb)
		no_annotation.write('\n')
		no_annotation.close
		record = SeqIO.read(gb, "genbank")
		gb_fasta = orf_dir + record.id + ".fasta" # 20160125
		with open(gb_fasta, 'w') as gbfasta:
			SeqIO.write(record, gbfasta, 'fasta')
		prodigal_command = ['prodigal', '-i', gb_fasta, '-d', orf_fasta, '-q']
		res = subprocess.call(prodigal_command)
		return
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

def draft_orf(strain, draft):
	print(draft)
	orf_fasta = orf_dir + strain + '_orf.fna' # change dir for orf fasta file 20180123
	prodigal_command = ['prodigal', '-i', draft, '-d', strain + '_temp.fna', '-q']
	res = subprocess.call(prodigal_command)
	with open(orf_fasta, 'w') as outfasta, open(strain + '_temp.fna', 'rU') as infasta:
		seq = SeqIO.parse(infasta, "fasta")
		for each_seq in seq:
			each_seq.id = strain + "_" + each_seq.id
			SeqIO.write(each_seq, outfasta, "fasta")
	os.remove(strain + '_temp.fna')

def orf_check(orf, newfile): #20180123
	blast_temp = tempfile.NamedTemporaryFile()
	db_name = blast_temp.name
	blast_temp.close
	mkbl_command = ['makeblastdb', '-dbtype', 'nucl', '-hash_index', \
	'-in', orf, '-out', db_name]
	res = subprocess.call(mkbl_command)
	blastn_cline = NcbiblastnCommandline(query=orf, db=db_name, evalue=ev,\
	outfmt=6, perc_identity=identity, qcov_hsp_perc=min_cover)
	stdout, stderr = blastn_cline()
	blast_results = stdout.rstrip().split('\n')
	stdout=''
	os.remove(db_name + '.nhr')
	os.remove(db_name + '.nin')
	os.remove(db_name + '.nog')
	os.remove(db_name + '.nsd')
	os.remove(db_name + '.nsi')
	os.remove(db_name + '.nsq')
	os.remove(db_name + '.nhd')
	os.remove(db_name + '.nhi')
	
	orf_tempdic = {}
	i = 0
	with open(orf, 'rU') as orf:
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
	with open(newfile, 'w') as newfasta:
		for each_seq, values in sorted(orf_tempdic.items(), key=lambda x: x[1][4]):
			if values[0] == 1:
				seq_items = SeqRecord(values[3], id = values[1], description = values[2])
				SeqIO.write(seq_items, newfasta, 'fasta')


def pickup_unique(strain_id):
	orf_tempdic = {}
	if seq_list[strain_id][2] == 'gb':
		gb_orf(seq_list[strain_id][1])
	elif seq_list[strain_id][2] == 'fasta':
		draft_orf(seq_list[strain_id][0], seq_list[strain_id][1])
	
	orf_fasta = orf_dir + seq_list[strain_id][0] + '_orf.fna' # change dir of orf fasta file 20180123
	nonredundant_fasta = orf_dir + seq_list[strain_id][0] + '_nonredundant.fna' # change dir of orf fasta file 20180123
	orf_check(orf_fasta, nonredundant_fasta)
	nonredundant_num = strain_id * 100000
	
	if strain_id == 0: #20180123
		with open(nonredundant_fasta, 'rU') as orf:
			sequence = SeqIO.parse(orf, 'fasta')
			for each in sequence:
				orf_dic[nonredundant_num] = [each.id, each.description, each.seq]
				nonredundant_num += 1
		return

	with open(nonredundant_fasta, 'rU') as orf:
		sequence = SeqIO.parse(orf, 'fasta')
		for each in sequence:
			orf_tempdic[each.id] = [0, each.id, each.description, each.seq]
	
	seqid_temp = tempfile.NamedTemporaryFile()
	with open(seqid_temp.name, 'w') as use_seq:
		for i in range(strain_id): # remove strain_id sequence from blastdb 180123
			use_seq.write(seq_list[i][0] + "\n")
		
	blastn_cline = NcbiblastnCommandline(query=nonredundant_fasta, db=out_name, evalue=ev,\
	outfmt=6, seqidlist=seqid_temp.name, perc_identity=identity, qcov_hsp_perc=min_cover,\
	num_alignments=2)
	stdout, stderr = blastn_cline()
	if stdout == '':
		for each_seq, values in orf_tempdic.items():
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

def make_binary(strain_id):
	print(seq_list[strain_id][0])
	seqid_temp = tempfile.NamedTemporaryFile()
	with open(seqid_temp.name, 'w') as use_seq:
			use_seq.write(seq_list[strain_id][0])
	blastn_cline = NcbiblastnCommandline(query=query_temp.name, db=out_name, evalue=ev,\
	outfmt=6, seqidlist=seqid_temp.name, perc_identity=identity, qcov_hsp_perc=min_cover,\
	num_alignments=1)
	stdout, stderr = blastn_cline()
	blast_results = stdout.rstrip().split('\n')
	stdout=''
	seqid_temp.close()
	for each_r in blast_results:
		br = each_r.split('\t')
		positive_list.append([strain_id, int(br[0])])		

def group_nex(fnum, glist, group_factor, join_factor):
	if fnum == 999:
		g_nex_filename = random_sample_nex
	else:
		g_nex_filename = out_name + '_group' + str(fnum).zfill(3) + '.nex' 
	temp_binary = binary_ndary[glist,:].astype(np.int)  ###### binary of selected strains
	orf_excluded_num = np.sum(temp_binary, axis = 0).astype(np.int)
	exclude_orf = []
	for i in range(orf_excluded_num.shape[0]):
		if orf_excluded_num[i] == 0:
			exclude_orf.append(i)
	temp_binary = np.delete(temp_binary, exclude_orf, 1)
	print('Writing binary sequence : ' + g_nex_filename)
	nexus(g_nex_filename, temp_binary, glist, 1, group_factor, join_factor)

def nexus(out_file, binaryseqs, seqlist, group_flag, group_factor, join_factor):
	
	strain = len(seqlist) #binaryseqs.shape[0]
	orf = binaryseqs.shape[1]
	div_lines = int(math.ceil(orf / 60.0))

	contents = open(out_file, 'w')
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
		contents.write(seq_list[seqlist[genome_id]][0] + "\n")
		for i in range(div_lines): #20180129 ##textwrap is very slow in cases of large sequences
			sequence = ""
			seq_start = i * 60
			seq_end = seq_start + 60
			if seq_end >= orf:
				seq_end = orf
			for each_binary in range(seq_start, seq_end):
				sequence = sequence + str(int(binaryseqs[genome_id, each_binary]))
			contents.write(sequence + "\n")
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
	for eachfile in seqlist:
		contents.write(seq_list[eachfile][0] + "\t" + seq_list[eachfile][1] + "\n")
	contents.write("------------------------------------------------------------\n")
	contents.write("]\n") # end comments

	
	contents.flush()
	contents.close

def calc_homology(seq1):
	second_seq = seq1 + 1
	dice_temp_list = [1 - distance.dice(binary_ndary[seq1],binary_ndary[i]) for i in range(num_strain)]
	dice_dic[seq1] = dice_temp_list
	share_temp_list = [np.sum(np.logical_and(binary_ndary[seq1], binary_ndary[i])) for i in range(num_strain)]
	share_dic[seq1] = share_temp_list

def part_share_indices(index_list):
	share_ndary = share_indices_ndary[index_list,:]
	share_ndary = share_ndary[:,index_list]
	return(share_ndary)

def grouping(sf, jf):
	remain_seq_list = range(num_strain)
	remain_seq_2nd = range(num_strain)
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

def out_mapping(group_list, out_fh):
	print('Making mapping file')
	temp_ndary = binary_ndary[group_list]
	temp_num_orf = temp_ndary.shape[1]
	temp_num_strain = temp_ndary.shape[0]
	del_bin_list = []
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

# br[0] query acc.ver
# br[1] subject acc.ver
# br[2] % identity
# br[3] alignment length
# br[4] mismatches
# br[5] gap opens
# br[6] q. start
# br[7] q. end
# br[8] s. start
# br[9] s. end
# br[10] evalue
# br[11] bit score

if __name__ == '__main__':
	if gbfiles == '' and fastafiles == '' and bin_file == '':
		print('Designate "gb file list" and/or "fasta file list", or "prebuild binary sequence csv.')
		sys.exit()

	start = time.time()
	out_fasta = out_name + "_elementORF.fasta"
	all_fasta = out_name + '_all.fasta'
	error_file = out_name + "_no_annotation.txt"
	current_dir = os.getcwd() # set current directory 20180123
	
	seq_dict = {}
	manager = Manager()
	seq_list = []
	acc_num_list =[]
	orf_dic = manager.dict() # unique orf contaner 20180124
	
	if bin_file != '':
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

	if gbfiles != '' and bin_file == '':
		file = open(gbfiles, "r")
		all_list = file.read().rstrip().split("\n")
		file.close
		print('Reading genbank files')
		while all_list.count("") >0:
			all_list.remove("")
		for each_gb in all_list:
			print(each_gb)
			record = SeqIO.read(each_gb, "genbank")
			acc_num = record.id.split('.')[0] # Jun. 13, 2018.
			if acc_num in acc_num_list:
				print(each_gb + ' is overlapped. Check and try again')
				sys.exit()
			if len(record.seq) < 20000:
				print(each_gb + ' is smaller than 20000. Remove and try again')
				os.remove(all_fasta)
				sys.exit()
			if len(record.seq) >= 1000000:
				grouping_flag = False
				
			acc_num_list.append(acc_num)
			seq_list.append([record.id, each_gb, 'gb'])
			with open(all_fasta, 'a') as longfasta:
				SeqIO.write(record, longfasta, 'fasta')
				longfasta.flush()

	if fastafiles != "" and bin_file == '':
		file = open(fastafiles, "r")
		all_draft = file.read().rstrip().split("\n")
		file.close
		print('Reading fasta files')
		while all_draft.count("") >0:
			all_draft.remove("")
		for each_draft in all_draft:
			draft_element = each_draft.split("\t")
			strain_name = draft_element[0]
			draft_path = draft_element[1]
			seq_list.append([strain_name, draft_path, 'fasta'])
			print(strain_name)
			draft_seqs = SeqIO.parse(draft_path, 'fasta')
			long_seq = ''
			for each_seq in draft_seqs:
				long_seq = long_seq + str(each_seq.seq) + \
				'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn'
			join_draft = SeqRecord(Seq(long_seq, IUPAC.ambiguous_dna), id = strain_name, description = '')
			with open(all_fasta, 'a') as longfasta:
				SeqIO.write(join_draft, longfasta, 'fasta')
				longfasta.flush()

	num_strain = len(seq_list)

	if bin_file == '':
		print('file summary writing')
	
		mkbl_command = ['makeblastdb', '-dbtype', 'nucl', '-parse_seqids', \
		'-in', all_fasta, '-out', out_name]
		res = subprocess.call(mkbl_command)
	
		orf_list = []
		query_temp = tempfile.NamedTemporaryFile()
		orf_path = out_name + '_orf'
		if os.path.isdir(orf_path): #20180123
			print('Create fasta files in "orf" directory')
		else:
			print('Create "orf" directory')
			os.mkdir(orf_path)
			print('Create fasta files in "orf" directory')
		orf_dir = current_dir + '/' + orf_path + '/'
		pool = Pool(num_process)	
		pool.map(pickup_unique, range(num_strain))
		pool.close()
		num_orf = len(orf_dic)
	
		print('Writing elemental ORFs') # 20180125
		with open(out_fasta, 'w') as outfh, open(query_temp.name, 'w') as queryfh:
			i = 0
			for key, values in sorted(orf_dic.items(), key=lambda x: x[0]):
				seq_items = SeqRecord(values[2], id = values[0], description = values[1])
				SeqIO.write(seq_items, outfh, 'fasta')
				query_items = SeqRecord(values[2], id = str(i))
				SeqIO.write(query_items, queryfh, 'fasta')
				orf_list.append(values[0])
				i += 1
		############# making binary sequences #################################
		print('Making binary sequences')
		positive_list = manager.list()
		pool = Pool(num_process)	
		pool.map(make_binary, range(num_strain))
		pool.close()
		query_temp.close()
		print('Constructiong binary array')
		binary_ndary = np.zeros((num_strain, num_orf), dtype = int) ##### create binary sequence array
		for values in positive_list:
			binary_ndary[values[0]][values[1]] = 1

		os.remove(out_name + '.nhr')
		os.remove(out_name + '.nin')
		os.remove(out_name + '.nog')
		os.remove(out_name + '.nsd')
		os.remove(out_name + '.nsi')
		os.remove(out_name + '.nsq')

	orf_dic.clear

	############# writing mapping file ####################################
	if mapping and bin_file == '':
		print('Making mapping file')
		num_orf = binary_ndary.shape[1]
		del_bin_list = []
		contents = open(out_name + '_mapping.txt', 'w')
		contents.write("\t")
		for genome_id in range(num_strain):
			contents.write(seq_list[genome_id][0] + "\t")
		contents.write("sum\n")
		for each_binary in range(num_orf):
			orf_sum = 0
			contents.write(orf_list[each_binary] + "\t")
			for genome_id in range(num_strain):
				contents.write(str(binary_ndary[genome_id, each_binary]) + "\t")
				orf_sum = orf_sum + binary_ndary[genome_id, each_binary]
			contents.write(str(orf_sum)+ '\n')
		contents.flush()
		contents.close
		print('Done! Mapping file : ' + out_name + '_mapping.txt')
	
	############# Calculating homology #############################
	print('Calculating homologies')
	share_dic = manager.dict()
	dice_dic = manager.dict()
	pool = Pool(num_process)	
	pool.map(calc_homology, range(num_strain))
	pool.close()
		
	share_indices_ndary = np.empty((0, num_strain), int)
	index_no = 0
	for i in range(num_strain):
		share_indices_ndary = np.append(share_indices_ndary, [share_dic[i]], axis = 0)

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
					if len(each_g) >= num_strain * 0.25:
						temp_share_index = part_share_indices(each_g)
						share_median_ndary = np.median(temp_share_index, axis = 0)
						num_lowsim = np.sum(share_median_ndary < 10)
						if num_lowsim > threshold:
							larger_num = max([num_lowsim, larger_num])
							print(larger_num)
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
					'_random' + str(random_sample) + '.nex'
					randomselect = random.sample(group_list[i], random_sample)
					group_nex(999, randomselect, g_factor, j_factor)#correct 181115
					with open(random_sample_nex, 'a') as outfh:
						outfh.write("[!\n")
						outfh.write("------- random sampling ------------------------------------\n")
						outfh.write(str(random_sample) + " sequences are selected randomly from group" + \
						str(i).zfill(3) + "\n")
						outfh.write("]\n")
				i += 1
 
	############# Writing whole binary sequence info #######################
	if bin_file == '':
		############# Output Nexus file ###########################
		seq_num_list = range(num_strain)
		print('Writing whole binary sequence')
		out_nexus = out_name + '_all.nex'
		nexus(out_nexus, binary_ndary, seq_num_list, 0, 0, 0)
		print('Done! Binary sequence file : ' + out_nexus)

		if bin_file == '':
			np.savetxt(out_name + '_binary.csv', binary_ndary, delimiter=',')
			seq_summary_file = out_name + '_info.csv'
			with open(seq_summary_file, 'w') as outfh:
				writer = csv.writer(outfh)
				for each in seq_list:
					writer.writerow(each)
			orf_list_file = out_name + '_orfinfo.txt'
			with open(orf_list_file, 'w') as outfh:
				for each in orf_list:
					outfh.write(each + '\n')

		############# Output Dice indices #########################
		dice_out = out_name + "_dice.txt"
		size = binary_ndary.shape
		seq_length = size[1]
		dice_line = ''
		for seq_num in range(len(seq_list)):
			dice_line = dice_line  + '\t'+ seq_list[seq_num][0]
		dice_list = [dice_line]
		d_indices = np.array([])
		for i in range(num_strain):
			mapped_dice = map(str, dice_dic[i])
			dice_line = '\t'.join(mapped_dice)
			dice_line = seq_list[i][0] + '\t' + dice_line
			dice_list.append(dice_line)
			d_indices = np.append(d_indices, np.array(dice_dic[i])) 

		with open(dice_out, 'w') as outfh:
			for line in dice_list:
				outfh.write(line + '\n')
			outfh.flush()
		print('Done! Dice indeces file : ' + dice_out)

		############# Output Share indices #########################
		share_out = out_name + "_share.txt"
		share_line = ''
		for seq_num in range(len(seq_list)):
			share_line = share_line  + '\t'+ seq_list[seq_num][0]
		share_list = [share_line]
		for i in range(num_strain):
			mapped_share = map(str, share_dic[i])
			share_line = '\t'.join(mapped_share)
			share_line = seq_list[i][0] + '\t' + share_line
			share_list.append(share_line)

		with open(share_out, 'w') as outfh:
			for line in share_list:
				outfh.write(line + '\n')
			outfh.flush()
		print('Done! Share indeces file : ' + share_out)

	############# Closing ##################################################
	elapsed_time = time.time() - start
	print("elapsed_time:{0}".format(elapsed_time) + "[sec]")
