# -*- coding: utf-8 -*-

"""
Created on Fri August 2 11:17:02 2024

@authors: Jon Pitsch and Martin Breuss
"""

import pandas as pd
import sys
import pysam
import os
import subprocess
from multiprocessing import Pool
import time
import math


dist = 95
extended = 250

def add_600bp(df, reference):
    '''get 300bp upstream and downstream of the position given in the df.
    Needs a pysam reference to work. fetch command appears to ignore the start
    number, but includes the end.'''
    
    seqs = []
    right_skwed = []
    left_skwed = []
    both_extended = []


    for index, row in df.iterrows():
        
        if row['CHROM'] in ['X', 'Y']:
            seq = reference.fetch('chr' + row['CHROM'], int(row['POS']) - dist - 1,
                                  int(row['POS']) + dist).upper()
            right = reference.fetch('chr' + row['CHROM'], int(row['POS']) - dist - 1,
                                  int(row['POS']) + dist + extended).upper()
            left = reference.fetch('chr' + row['CHROM'], int(row['POS']) - dist - extended - 1,
                                  int(row['POS']) + dist).upper()
        else:		
            seq = reference.fetch('chr' + str(int(row['CHROM'])), int(row['POS']) - dist - 1,
                                  int(row['POS']) + dist).upper()
            right = reference.fetch('chr' + str(int(row['CHROM'])), int(row['POS']) - dist - 1,
                                  int(row['POS']) + dist + extended).upper()
            left = reference.fetch('chr' + str(int(row['CHROM'])), int(row['POS']) - dist - extended - 1,
                                  int(row['POS']) + dist).upper()

        seqs.append(seq)
        right_skwed.append(right)
        left_skwed.append(left)
        
    df['Sequence_for_Primer3'] = pd.Series(seqs, index=df.index)
    df['Sequence_for_Primer3_right'] = pd.Series(right_skwed, index=df.index)
    df['Sequence_for_Primer3_left'] = pd.Series(left_skwed, index=df.index)

standard_primer_names = []
right_primers = []
left_primers = []
all_primer_names = []
on_target_pos = []

#def write_primer3_file(df, path_to_output, padding=60):
def write_primer3_file(df, path_to_output, padding):
    '''using the df with the 600bp it will write a primer3 file that will allow
    to batch design suitable primer pairs. padding defines the upstream and
    downstream padding on each site that should be used. path_to_thermo has to
    lead to the folder that contains the primer3_configstack.ds file. trailing
    / is added by program.'''
    
    start = dist - padding
    start_left = dist + extended - padding
    pad2 = padding * 2
    extended_padding = padding * 3
    
    with open(path_to_output, 'a') as output:
        
        for index, row in df.iterrows():

            output.write(('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n' +
                          #'PRIMER_TASK=generic\n' +
                          #'PRIMER_PICK_LEFT_PRIMER=1\n' +
                          #'PRIMER_PICK_INTERNAL_OLIGO=0\n' +
                          #'PRIMER_PICK_RIGHT_PRIMER=1\n' +
                          #'PRIMER_OPT_SIZE=20\n' +
                          #'PRIMER_MIN_SIZE=18\n' +
                          #'PRIMER_MAX_SIZE=22\n' +
                          'SEQUENCE_TARGET={},{}\n' +
                          'PRIMER_PRODUCT_SIZE_RANGE=100-500\n' +
                          'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=./primer3_config\n' +
                          #'PRIMER_EXPLAIN_FLAG=1\n' +
                          '=\n').format(str(row['PROJ']) + '_' + str(row['CHROM']) + '_' +
                                       str(row['POS']) + '_ideal',
                                       row['Sequence_for_Primer3'], str(start),
                                       str(pad2)))

            standard_primer = str(row['PROJ']) + '_' + str(row['CHROM']) + '_' + str(row['POS']) + '_ideal' 

            standard_primer_names.append(standard_primer)

            all_primer_names.append(standard_primer)

def write_primer3_file_right(df, path_to_output, padding):
    '''using the df with the 600bp it will write a primer3 file that will allow
    to batch design suitable primer pairs. padding defines the upstream and
    downstream padding on each site that should be used. path_to_thermo has to
    lead to the folder that contains the primer3_configstack.ds file. trailing
    / is added by program.'''
    
    start = dist - padding
    start_left = dist + extended - padding
    pad2 = padding * 2
    extended_padding = padding * 3
    
    with open(path_to_output, 'a') as output:
        
        for index, row in df.iterrows():

            output.write(('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n' +
                          #'PRIMER_TASK=generic\n' +
                          #'PRIMER_PICK_LEFT_PRIMER=1\n' +
                          #'PRIMER_PICK_INTERNAL_OLIGO=0\n' +
                          #'PRIMER_PICK_RIGHT_PRIMER=1\n' +
                          #'PRIMER_OPT_SIZE=20\n' +
                          #'PRIMER_MIN_SIZE=18\n' +
                          #'PRIMER_MAX_SIZE=22\n' +
                          'SEQUENCE_TARGET={},{}\n' +
                          'PRIMER_PRODUCT_SIZE_RANGE=100-500\n' +
                          'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=./primer3_config\n' +
                          #'PRIMER_EXPLAIN_FLAG=1\n' +
                          '=\n').format(str(row['PROJ']) + '_' + str(row['CHROM']) + '_' +
                                       str(row['POS']) + '_right',
                                       row['Sequence_for_Primer3_right'], str(start),
                                       str(extended_padding)))
            
            right_primer = str(row['PROJ']) + '_' + str(row['CHROM']) + '_' + str(row['POS']) + '_right'

            right_primers.append(right_primer)

            all_primer_names.append(right_primer)

def write_primer3_file_left(df, path_to_output, padding):
    '''using the df with the 600bp it will write a primer3 file that will allow
    to batch design suitable primer pairs. padding defines the upstream and
    downstream padding on each site that should be used. path_to_thermo has to
    lead to the folder that contains the primer3_configstack.ds file. trailing
    / is added by program.'''
    
    start = dist - padding
    start_left = dist + extended - padding
    pad2 = padding * 2
    extended_padding = padding * 3
    
    with open(path_to_output, 'a') as output:
        
        for index, row in df.iterrows():
        
            output.write(('SEQUENCE_ID={}\nSEQUENCE_TEMPLATE={}\n' +
                          #'PRIMER_TASK=generic\n' +
                          #'PRIMER_PICK_LEFT_PRIMER=1\n' +
                          #'PRIMER_PICK_INTERNAL_OLIGO=0\n' +
                          #'PRIMER_PICK_RIGHT_PRIMER=1\n' +
                          #'PRIMER_OPT_SIZE=20\n' +
                          #'PRIMER_MIN_SIZE=18\n' +
                          #'PRIMER_MAX_SIZE=22\n' +
                          'SEQUENCE_TARGET={},{}\n' +
                          'PRIMER_PRODUCT_SIZE_RANGE=100-500\n' +
                          'PRIMER_THERMODYNAMIC_PARAMETERS_PATH=./primer3_config\n' +
                          #'PRIMER_EXPLAIN_FLAG=1\n' +
                          '=\n').format(str(row['PROJ']) + '_' + str(row['CHROM']) + '_' +
                                       str(row['POS']) + '_left',
                                       row['Sequence_for_Primer3_left'], str(start_left),
                                       str(padding)))   
            
            left_primer = str(row['PROJ']) + '_' + str(row['CHROM']) + '_' + str(row['POS']) + '_left'

            left_primers.append(left_primer)

            all_primer_names.append(left_primer)

def parse_p3_output(filepath, counterflag=True, pairs=2):
    '''Taking the machine output of primer3 command line tool, function parses
    the file for the information line and the lines containing the primer
    sequences. If the counterflag is True, then it will stop after obtaining
    the number of pairs equaling pairs (2 by default). Deactivating pairs will
    return all 5 primers that are returned as standard output. '=' will reset
    counter, so this should also work if there is only less pairs than defined
    by pairs.'''
    names = []
    name_worked = []
    
    forward_sequences = []
    reverse_sequences = []

    forward_tm = []
    reverse_tm = []

    print('Reading P3 Output')

    with open(filepath) as file:
        name = ''
        counter = 0

        for line in file:

            line_lst = [element.split('_') for element in line.split('=')]

            if line == '=\n':
                counter = 0

            elif counter == (pairs*2) and counterflag == True:
                continue

            elif line.startswith('SEQUENCE_ID'):
                name = line.split('=')[1].replace('\n', '')


            elif len(line_lst[0]) > 3 and counter == 3:
                if line_lst[0][1] == 'RIGHT' and line_lst[0][3] == 'TM':
                    # Append twice. Once for forward-forward and once for actual primer name
                    reverse_tm.append(line_lst[1][0].replace('\n', ''))
                    forward_tm.append(line_lst[1][0].replace('\n', ''))
                    reverse_tm.append(line_lst[1][0].replace('\n', ''))
                    counter += 1
                    continue

            elif len(line_lst[0]) > 3 and counter == 2:
                if line_lst[0][1] == 'LEFT' and line_lst[0][3] == 'TM':
                    # Append twice. Once for forward-forward and once for actual primer name
                    forward_tm.append(line_lst[1][0].replace('\n', ''))
                    reverse_tm.append(line_lst[1][0].replace('\n', ''))
                    forward_tm.append(line_lst[1][0].replace('\n', ''))
                    counter += 1
                    continue
                else:
                    if not_present==1:
                        tm = '0.0'
                        reverse_tm.append(tm)
                        forward_tm.append(tm)
                        counter += 1
                        continue

            elif len(line_lst[0]) > 3 and counter == 1:
                if line_lst[0][3] == 'SEQUENCE':
                    reverse_sequences.append(line_lst[1][0].replace('\n', ''))
                    names.append('-'.join([name,'reverse-reverse']))
                    name_worked.append('-'.join([name,'reverse-reverse']))
                    forward_sequences.append(line_lst[1][0].replace('\n', ''))
                    reverse_sequences.append(line_lst[1][0].replace('\n', ''))
                    counter += 1
                    continue

            elif len(line_lst[0]) > 3 and counter == 0:
                if line_lst[0][3] == 'SEQUENCE':
                    names.append('-'.join([name,'forward-forward']))
                    name_worked.append('-'.join([name,'forward-forward']))
                    forward_sequences.append(line_lst[1][0].replace('\n', ''))
                    reverse_sequences.append(line_lst[1][0].replace('\n', ''))
                    names.append(name)
                    name_worked.append(name)
                    forward_sequences.append(line_lst[1][0].replace('\n', ''))
                    counter += 1
                    continue
                elif line_lst[0][1] == 'PAIR' and line_lst[0][3] == 'RETURNED':
                    line_lst[1][0].replace('\n', '')
                    line_lst[1][0] = int(line_lst[1][0])
                    if line_lst[1][0] == 0:
                        not_present = 1
                        names.append(name)
                        seq = 'AAAAAAAAAAAAAAAAAAA'
                        forward_sequences.append(seq)
                        reverse_sequences.append(seq)
                        counter += 1
                        continue
                else:
                    continue
            else:
                continue

    print(len(name_worked))
    print(len(forward_tm))

    df = pd.DataFrame(names, columns=['Primer_Name'])
    df['Forward_Primer'] = pd.Series(forward_sequences, index=df.index,dtype='object')
    df['Reverse_Primer'] = pd.Series(reverse_sequences, index=df.index,dtype='object')

    df2 = pd.DataFrame(name_worked, columns=['Primer_Name'])
    df2['Forward_Tm'] = pd.Series(forward_tm, index=df2.index,dtype='float')
    df2['Reverse_Tm'] = pd.Series(reverse_tm, index=df2.index,dtype='float')

    df_merged = df.merge(df2,how='left',on='Primer_Name')
    df_merged['Forward_Tm'] = df_merged['Forward_Tm'].fillna(0)
    df_merged['Reverse_Tm'] = df_merged['Reverse_Tm'].fillna(0)

    return df_merged

# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# execution

print(
'''Usage of {0}: python3 {0} data_path ref_path therm_path [padding [pairs]

input:
    
data: .csv file; columns labeled 'CHROM' and 'POS'
ref: genome reference used for positional information in data; .fasta file
therm: path to folder that contains the primer3_configstack.ds

optional:
padding: int; # of bp up- and downstream of variant position; default=50
pairs: int; # of primer pairs returned; default=2, max=5 (primer3 setting)

if pairs is defined, padding has to be defined too!


output (in current folder):
    
csv with 601bp sequence
primer3 input file
human readable primer3 output
machine readable primer3 output
csv with primer names and sequences.
'''.format(sys.argv[0]))

#cont = input('Do you want to continue? [y/n] ')

#if cont.upper() not in ['YES', 'Y']:  
#    sys.exit()
    
#~~~~~
#input check
    
if len(sys.argv) < 5 or len(sys.argv) > 8:
    print('Wrong number of arguments! See usage above.')
    sys.exit()

if os.path.splitext(sys.argv[1])[1] != '.csv':
    print('Requires .csv as data input.')
    sys.exit()

#input check done
#~~~~~~
    
data_path = sys.argv[1]
ref_path = sys.argv[2]
reference = pysam.Fastafile(ref_path)
#therm_path = sys.argv[3].rstrip('/')
output_dir = sys.argv[3]
name = sys.argv[4]


if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not output_dir.endswith("/"):
    output_dir += "/"

global padding 
padding = 60
pairs = 2
print(len(sys.argv))
#if len(sys.argv) > 7:
#    padding = int(sys.argv[7])
#    print(padding)
#if len(sys.argv) == 8:
#    pairs = int(sys.argv[8])

data = pd.read_csv(data_path)

if int(len(data['CHROM'])) > 100000:
    batch_size = int(len(data['CHROM'])/10)
else:
    batch_size = int(len(data['CHROM']))


###### Primer3 Steps for Ideal ####
primer_start_time = time.time()

add_600bp(data, reference)

type_string = '_ideal'

data.to_csv(name + '_601bp_region'+type_string+'.csv', index=False)

write_primer3_file(data, name + '_primer3_input'+type_string+'.txt',
                   60)

with open(name + '_primer3_input'+type_string+'.txt', 'rb', 0) as input_f, \
open(name + '_human_output'+type_string+'.txt', 'wb', 0) as output_f:
    subprocess.run(['primer3_core', '-format_output'], stdin=input_f,
                   stdout=output_f, check=True)

with open(name + '_primer3_input'+type_string+'.txt', 'rb', 0) as input_f, \
open(name + '_machine_output'+type_string+'.txt', 'wb', 0) as output_f:
    subprocess.run(['primer3_core'], stdin=input_f,
                   stdout=output_f, check=True)

output = parse_p3_output(name + '_machine_output'+type_string+'.txt', pairs=pairs)

output.to_csv(name + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

print("Ideal Primers Created in --- %s seconds ---" % (time.time() - primer_start_time))
#### End of Primer3 Steps for Ideal #####

#all_primers = pd.DataFrame()
#all_primers['primers'] = all_primer_names
#all_primers.to_csv('all_primers_list.csv',index=False)

row_split = []

def round_up(n, decimals = 0): 
    multiplier = 10 ** decimals 
    return math.ceil(n * multiplier) / multiplier

def split_primer_list():
    #row_split = []
    primers = pd.read_csv(''+name+'_primer_list'+type_string+'.txt', delimiter='\t',names = ['primer','forward','reverse', 'forward_tm', 'reverse_tm'])

    x = len(primers['primer'])

    y = int(round_up(x, -2))

    split = []

    for i in range(0,y+batch_size,batch_size):
        iteration = i
        split.append(iteration)

    if x < batch_size:
        row_split.append(x)
        print(row_split)
        primers.to_csv(name+'_'+str(x) + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

    print('This is split')
    print(split)

    if x > batch_size:
        for i in split:
            if i == 0:
                continue
            if i == batch_size:
                row = i
                row_split.append(row)
                df = primers.iloc[:batch_size,:]
                df.to_csv(name+'_'+str(i) + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

            elif i > batch_size and x > i:
                row = i
                row_split.append(row)
                df = primers.iloc[i-batch_size:i]
                df.to_csv(name+'_'+str(i) + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

            elif i > batch_size and x < i and x-i > -1*(batch_size-1):
                row = x
                row_split.append(row)
                df = primers.iloc[i-batch_size:,:]
                df.to_csv(name+'_'+str(x) + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

            else:
                break

split_primer_list()

print(row_split)

def ispcr():
    #print(row_split)
    #print(type_string)

    print('Running isPcr')

    #start_time = time.time()

    for i in row_split:

        start_time = time.time()

        subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list'+type_string+'.txt stdout -out=bed -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr'+type_string+'.bed',shell=True)

        subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list'+type_string+'.txt stdout -out=fa -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr'+type_string+'.fa',shell=True)

        print("--- %s seconds ---" % (time.time() - start_time))

ispcr()
#pool = Pool(8)
#res= pool.map(ispcr, row_split, 8)

ideal_total_time = time.time()

def ispcr_analysis():

    total_start_time = time.time()
    print('Running P3-isPcr analysis')

    for f in row_split:
        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:

            continue
        ispcr = pd.read_csv(name+'_'+str(f)+'_isPcr'+type_string+'.bed',delimiter='\t',names=['chrom','amplicon_start','amplicon_end','primer_name','score','strand'])

        def fasta_reader(filename):
            from Bio.SeqIO.FastaIO import FastaIterator
            with open(filename) as handle:
                for record in FastaIterator(handle):
                    yield record

        sequences = []

        for entry in fasta_reader(""+name+"_"+str(f)+"_isPcr"+type_string+".fa"):
            seq = str(entry.seq)
            sequences.append(seq)

        chrom = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

        ispcr['amplicon'] = sequences

        chromo = []
        amp_start = []
        amp_end = []
        primer = []
        score = []
        size = []
        sequence = []
        ID = []

        for i in range(len(ispcr['chrom'])):
            if ispcr['chrom'][i] in chrom:
                chro = ispcr['chrom'][i]
                start = ispcr['amplicon_start'][i]
                end = ispcr['amplicon_end'][i]
                sz = abs((ispcr['amplicon_start'][i])-(ispcr['amplicon_end'][i]))
                prim = ispcr['primer_name'][i]
                sco = ispcr['score'][i]
                seq = ispcr['amplicon'][i]
                variant_id = str(chro) + '_' + str(start) + '_' + str(end) + '_' + str(sz) + '_' + str(prim) + '_' + str(sco)
            else:
                continue

            chromo.append(chro)
            amp_start.append(start)
            amp_end.append(end)
            primer.append(prim)
            score.append(sco)
            size.append(sz)
            sequence.append(seq)
            ID.append(variant_id)

        results = pd.DataFrame()
        results['chrom'] = chromo
        results['amplicon_start'] = amp_start
        results['amplicon_end'] = amp_end
        results['amplicon_size'] = size
        results['primer_name'] = primer
        results['variant_ID'] = ID

        p3_output = pd.read_csv(name+'_'+str(f)+'_primer_list'+type_string+'.txt',delimiter='\t',names=['primer_name','forward_primer','reverse_primer','forward_tm', 'reverse_tm'])

        forward = []
        reverse = []

        forward_tm = []
        reverse_tm = []

        forw=0
        rev=0
        forw_tm = 0
        rev_tm = 0

        for k in range(len(primer)):
            for j in range(len(p3_output)):
                if p3_output['primer_name'][j] ==  primer[k]:
                    forw = p3_output['forward_primer'][j]
                    forw_tm = p3_output['forward_tm'][j]
                    rev = p3_output['reverse_primer'][j]
                    rev_tm = p3_output['reverse_tm'][j]
                else:
                    continue

            forward.append(forw)
            forward_tm.append(forw_tm)
            reverse.append(rev)
            reverse_tm.append(rev_tm)

        results['forward_name'] = results['primer_name'] + '_F'
        results['forward_primer'] = forward
        results['forward_tm'] = forward_tm
        results['reverse_name'] = results['primer_name'] + '_R'
        results['reverse_primer'] = reverse
        results['reverse_tm'] = reverse_tm
        results['score'] = score
    
        results['amplicon'] = sequence

        print("Finished ISPCR analysis in --- %s seconds ---" % (time.time() - total_start_time))

        results.to_csv(name+'_'+str(f)+'_p3_ispcr_output'+type_string+'.txt',sep='\t',index=False)

ispcr_analysis()

def filter_ispcr_analysis():

    total_start_time = time.time()
    print('Filtering isPcr results. Removing any degenerate primer with a score less than 750. Removing duplicate sites.')

    for i in row_split:
        print(i)

        if int(os.path.getsize(name+'_'+str(i)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        df = pd.read_csv(name+'_'+str(i)+'_p3_ispcr_output'+type_string+'.txt',delimiter='\t')

        df2 = df.drop_duplicates(subset=['variant_ID'])

        df2['good_score'] = (df2['score'] >= 750)

        df_true=df2[df2.pop('good_score')]
        
        df_true['on_target_chrom'] = 'chr'+df_true['primer_name'].str.split('_', expand=False).str[-3]
        df_true['on_target_pos'] = df_true['primer_name'].str.split('_', expand=False).str[-2]
        df_true['on_target_pos'] = df_true['on_target_pos'].astype(int)

        df_true.to_csv(name+'_'+str(i)+'_score_filtered'+type_string+'.txt',index=False,sep='\t')

        print("Finished Filtering ISPCR results in --- %s seconds ---" % (time.time() - total_start_time))

filter_ispcr_analysis()

def split_fasta():

    total_start_time = time.time()
    print('Splitting primers from amplicon sequence. Counting mismatches.')

    for f in row_split:
        
        forw_mis = []
        rev_mis = []
        lower_forward_count = []
        lower_reverse_count = []
        forward_mis_string = []
        reverse_mis_string = []

        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        df = pd.read_csv(name+'_'+str(f)+'_score_filtered'+type_string+'.txt',delimiter='\t')

        for i in range(len(df['amplicon'])):
            forw = len(df['forward_primer'][i])
            rev = len(df['reverse_primer'][i])

            reverse = len(df['amplicon'][i]) - rev

            forw_p, rev_p = df['amplicon'][i][:forw], df['amplicon'][i][reverse:]

            forw_mis.append(forw_p)
            rev_mis.append(rev_p)
        
            def n_lower_chars_forward(string):
                lower_f = sum(map(str.islower, string))
            
                lower_forward_count.append(lower_f)

            n_lower_chars_forward(forw_p)
        
            def n_lower_chars_reverse(string):
                lower_r = sum(map(str.islower, string))

                lower_reverse_count.append(lower_r)

            n_lower_chars_reverse(rev_p)

            rever_mis = '3\' '+str.join('', ('X' if chr.islower() else '_' for chr in rev_p))+' 5\''
        
            reverse_mis_string.append(rever_mis)

            forwa_mis = '5\' '+str.join('', ('X' if chr.islower() else '_' for chr in forw_p))+' 3\''

            forward_mis_string.append(forwa_mis)

        df['forward_alignment'] = forw_mis
        df['mismatch_count_forward'] = lower_forward_count
        df['forward_alignment_string'] = forward_mis_string
        df['reverse_alignment'] = rev_mis
        df['mismatch_count_reverse'] = lower_reverse_count
        df['reverse_alignment_string'] = reverse_mis_string

        df.to_csv(name+'_'+str(f)+'_final_output'+type_string+'.txt',sep='\t',index=False)
        
        print("Finished Finding Mismatches in --- %s seconds ---" % (time.time() - total_start_time))
        
split_fasta()

def standard_present():

    total_start_time = time.time()
    for f in row_split:

        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        df = pd.read_csv(name+'_'+str(f)+'_final_output'+type_string+'.txt',delimiter='\t')

        standard = []
        for i in range(len(df['primer_name'])):
            if df['primer_name'][i] in standard_primer_names:
                p_name = 'TRUE'
                standard.append(p_name)

            else:
                p_name = 'FALSE'
                standard.append(p_name)
  
        df['standard'] = standard

        df.to_csv(name+'_'+str(f)+'_standard_check'+type_string+'.txt',sep='\t',index=False)

        print("Finished Checking Standard Presence in --- %s seconds ---" % (time.time() - total_start_time))

standard_present()

from difflib import SequenceMatcher
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.gap_score = -0.10

def standard_matches():

    total_start_time = time.time()

    for f in row_split:

        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        df = pd.read_csv(name+'_'+str(f)+'_standard_check'+type_string+'.txt',delimiter='\t')
        good_match = []
        amplicons = []
        percent = []
        match_good_match = []
        align_percent = []
        align_percent_standard = []

        for i in range(len(all_primer_names)):
            for p in range(len(df['primer_name'])):
                if (df['primer_name'][p] == all_primer_names[i]) and (df['score'][p] == 1000) and (df['on_target_chrom'][p] == df['chrom'][p]) and (df['amplicon_start'][p] < df['on_target_pos'][p]) and (df['amplicon_end'][p] > df['on_target_pos'][p]):
                    good = df['primer_name'][p]
                    amp = df['amplicon'][p]
                    amplicons.append(amp)
                    good_match.append(good)

        df3 = pd.DataFrame()

        df3['good_match'] = good_match
        df3['amplicons'] = amplicons

        df2 = df3.drop_duplicates(subset=['good_match'])

        df2.to_csv(name+'_'+str(f)+'_good_match_and_amplicons'+type_string+'.txt', sep='\t', index=False)

        print("Finished Checking Standard Matches in --- %s seconds ---" % (time.time() - total_start_time))

standard_matches()

def percent_match():

    total_start_time = time.time()

    aligner = Align.PairwiseAligner()
    aligner.gap_score = -0.10

    for f in row_split:

        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        df = pd.read_csv(name+'_'+str(f)+'_standard_check'+type_string+'.txt',delimiter='\t')

        df2 = pd.read_csv(name+'_'+str(f)+'_good_match_and_amplicons'+type_string+'.txt',delimiter='\t')

        scores = [-1] * len(df['primer_name'])
        stand_scores = [-1] * len(df['primer_name'])

        df['percent_match'] = 'On_Target_Not_Found_by_isPcr'

        df['gold_amplicon'] = 'On_Target_Not_Found_by_isPcr'

        for j in range(len(df['primer_name'])):
            for k in range(len(df2['good_match'])):
                if df['primer_name'][j] == df2['good_match'][k]+'-forward-forward':
                    test_amp = df['amplicon'][j]
                    standard_amp = df2['amplicons'][k]

                    percent_match = SequenceMatcher(None, standard_amp, test_amp).ratio()

                    df['percent_match'][j] = percent_match
                
                    align = aligner.align(standard_amp, test_amp)

                    score = align.score/len(test_amp)

                    scores[j] = score

                    stand_score = align.score/len(standard_amp)

                    stand_scores[j] = stand_score

                    match = df2['good_match'][k]

                    df['gold_amplicon'][j] = match
                    continue
                elif df['primer_name'][j] == df2['good_match'][k]+'-reverse-reverse':
                    test_amp = df['amplicon'][j]
                    standard_amp = df2['amplicons'][k]

                    percent_match = SequenceMatcher(None, standard_amp, test_amp).ratio()
                    
                    df['percent_match'][j] = percent_match

                    align = aligner.align(standard_amp, test_amp)

                    score = align.score/len(test_amp)

                    scores[j] = score

                    stand_score = align.score/len(standard_amp)

                    stand_scores[j] = stand_score

                    match = df2['good_match'][k]

                    df['gold_amplicon'][j] = match
                    continue
                elif df['primer_name'][j] == df2['good_match'][k]:
                    test_amp = df['amplicon'][j]
                    standard_amp = df2['amplicons'][k]

                    percent_match = SequenceMatcher(None, standard_amp, test_amp).ratio()

                    df['percent_match'][j] = percent_match

                    align = aligner.align(standard_amp, test_amp)

                    score = align.score/len(test_amp)

                    scores[j] = score

                    stand_score = align.score/len(standard_amp)

                    stand_scores[j] = stand_score

                    match = df2['good_match'][k]

                    df['gold_amplicon'][j] = match
                    continue

        df['normalized_match_to_test_amplicon'] = scores

        df['normalized_match_to_gold_amplicon'] = stand_scores

        df.to_csv(name+'_'+str(f)+'_final_output_with_percent_match'+type_string+'.txt',sep='\t',index=False)

        print("Finished Calculating Percent Matches in --- %s seconds ---" % (time.time() - total_start_time))

percent_match()

input_list = pd.read_csv(data_path)
input_list['variant_id'] = input_list['PROJ'].astype(str) + '_' + input_list['CHROM'].astype(str) + '_' + input_list['POS'].astype(str)

primer_list = pd.read_csv(name+'_primer_list'+type_string+'.txt',delimiter='\t', names = ['primer_name', 'forward','reverse'])

ideal_primer = []
alternate_primer = []
forward_name = []
forward_primer_list = []
forward_tm_list = []
reverse_name = []
reverse_primer_list = []
reverse_tm_list = []
amplicon_name = []
off_target = []
amplicon_start = []
amplicon_end = []
amplicon_size = []
primer_count_high = []

def capture_best_matches():
    for f in row_split:

        if int(os.path.getsize(name+'_'+str(f)+'_isPcr'+type_string+'.bed')) == 0:
            continue

        p3i_output = pd.read_csv(name+'_'+str(f)+'_final_output_with_percent_match'+type_string+'.txt',delimiter='\t')

        primer_counts = []

        for j in range(len(p3i_output['gold_amplicon'])):
            count = len(p3i_output[p3i_output['gold_amplicon'] == p3i_output['gold_amplicon'][j]])
            primer_counts.append(count)
        
        p3i_output['primer_count'] = primer_counts

        def primer_evaluation(df):
            df['ideal_primer'] = ((df['score'].astype(int) ==1000) & (df['standard'].astype(str) == 'True')
                  & (df['normalized_match_to_test_amplicon'] == 1.000)
                  & (df['normalized_match_to_gold_amplicon'] == 1.000) & (df['on_target_chrom'] == df['chrom']) & (df['amplicon_start'] < df['on_target_pos']) & (df['amplicon_end'] > df['on_target_pos']))
        
            df['alternate_primer'] = ((df['score'].astype(int) ==1000) & (df['standard'].astype(str) == 'False')
                  & (df['normalized_match_to_test_amplicon'] == 1.000)
                  & (df['normalized_match_to_gold_amplicon'] == 1.000) & (df['on_target_chrom'] == df['chrom']) & (df['amplicon_start'] < df['on_target_pos']) & (df['amplicon_end'] > df['on_target_pos']))


        primer_list = pd.read_csv(name+'_'+str(f)+'_primer_list'+type_string+'.txt',delimiter='\t', names = ['primer_name', 'forward','reverse'])


        primer_evaluation(p3i_output)

        def off_target_check(df):
            df['high_con_off_target_1'] = ((df['normalized_match_to_test_amplicon'].astype(float) < 1.00) & (df['normalized_match_to_test_amplicon'].astype(float) > 0.80))

            df['high_con_off_target_2'] = ((df['normalized_match_to_gold_amplicon'].astype(float) < 1.00) & (df['normalized_match_to_gold_amplicon'].astype(float) > 0.80))

            df['perfect_offtarget'] = ((df['score'].astype(int) ==1000) & (df['ideal_primer'].astype(str) == 'False') & (df['alternate_primer'].astype(str) == 'False'))

        off_target_check(p3i_output)
        
        poss_off_targets = []

        for j in range(len(p3i_output['gold_amplicon'])):
            if p3i_output['primer_count'][j] >= 10:
                poss = 'True'
                poss_off_targets.append(poss)
                continue
            elif p3i_output['primer_count'][j] < 10:
                temp = p3i_output[p3i_output['gold_amplicon'] == p3i_output['gold_amplicon'][j]]
                high_con_1 = temp['high_con_off_target_1'].to_list()
                high_con_2 = temp['high_con_off_target_2'].to_list()
                perfect_match = temp['perfect_offtarget'].to_list()
                poss = 0
                for entry in range(len(high_con_1)):
                    if high_con_1[entry] == True:
                        poss = 'True'
                        poss_off_targets.append(poss)
                        break
                    elif high_con_2[entry] == True:
                        poss = 'True'
                        poss_off_targets.append(poss)
                        break
                    elif perfect_match[entry] == True:
                        poss = 'True'
                        poss_off_targets.append(poss)
                        break
                    else:
                        continue
                if poss == 0:
                    poss = 'False'
                    poss_off_targets.append(poss)
                    continue
            else:
                poss='False'
                poss_off_targets.append(poss)
                continue

        p3i_output['possible_off_targets'] = poss_off_targets

        for j in range(len(p3i_output['gold_amplicon'])):
            if (p3i_output['ideal_primer'][j].astype(str) == 'True'):
                ideal = 'True'
                alt = 'False'
                forward_n = p3i_output['forward_name'][j]
                forward_p = p3i_output['forward_primer'][j]
                forward_t = p3i_output['forward_tm'][j]
                reverse_n = p3i_output['reverse_name'][j]
                reverse_p = p3i_output['reverse_primer'][j]
                reverse_t = p3i_output['reverse_tm'][j]
                amp = p3i_output['gold_amplicon'][j]
                poss = p3i_output['possible_off_targets'][j]
                amp_start = p3i_output['amplicon_start'][j]
                amp_end = p3i_output['amplicon_end'][j]
                amp_size = p3i_output['amplicon_size'][j]
                count = p3i_output['primer_count'][j]
                ideal_primer.append(ideal)
                alternate_primer.append(alt)
                forward_name.append(forward_n)
                forward_primer_list.append(forward_p)
                forward_tm_list.append(forward_t)
                reverse_name.append(reverse_n)
                reverse_primer_list.append(reverse_p)
                reverse_tm_list.append(reverse_t)
                amplicon_name.append(amp)
                off_target.append(poss)
                amplicon_start.append(amp_start)
                amplicon_end.append(amp_end)
                amplicon_size.append(amp_size)
                primer_count_high.append(count)
                continue
            if (p3i_output['alternate_primer'][j].astype(str) == 'True'):
                ideal_alt = 'False'
                alt_alt = 'True'
                forward_n_alt = p3i_output['forward_name'][j]
                forward_p_alt = p3i_output['forward_primer'][j]
                forward_t_alt = p3i_output['forward_tm'][j]
                reverse_n_alt = p3i_output['reverse_name'][j]
                reverse_p_alt = p3i_output['reverse_primer'][j]
                reverse_t_alt = p3i_output['reverse_tm'][j]
                amp_alt = p3i_output['gold_amplicon'][j]
                poss_alt = p3i_output['possible_off_targets'][j]
                amp_start_alt = p3i_output['amplicon_start'][j]
                amp_end_alt = p3i_output['amplicon_end'][j]
                amp_size_alt = p3i_output['amplicon_size'][j]
                count_alt = p3i_output['primer_count'][j]
                ideal_primer.append(ideal_alt)
                alternate_primer.append(alt_alt)
                forward_name.append(forward_n_alt)
                forward_primer_list.append(forward_p_alt)
                forward_tm_list.append(forward_t_alt)
                reverse_name.append(reverse_n_alt)
                reverse_primer_list.append(reverse_p_alt)
                reverse_tm_list.append(reverse_t_alt)
                amplicon_name.append(amp_alt)
                off_target.append(poss_alt)
                amplicon_start.append(amp_start_alt)
                amplicon_end.append(amp_end_alt)
                amplicon_size.append(amp_size_alt)
                primer_count_high.append(count_alt)
                continue
            elif (p3i_output['ideal_primer'][j].astype(str) == 'False') & (p3i_output['alternate_primer'][j].astype(str) == 'False') & (p3i_output['on_target_chrom'][j] == p3i_output['chrom'][j]) & (p3i_output['amplicon_start'][j] < p3i_output['on_target_pos'][j]) & (p3i_output['amplicon_end'][j] > p3i_output['on_target_pos'][j]) & (p3i_output['score'][j] == 1000) & (p3i_output['normalized_match_to_test_amplicon'][j] > .97) & (p3i_output['normalized_match_to_gold_amplicon'][j] > .97):
                ideal_uniq = 'False'
                alt_uniq = 'False'
                forward_n_uniq = p3i_output['forward_name'][j]
                forward_p_uniq = p3i_output['forward_primer'][j]
                forward_t_uniq = p3i_output['forward_tm'][j]
                reverse_n_uniq = p3i_output['reverse_name'][j]
                reverse_p_uniq = p3i_output['reverse_primer'][j]
                reverse_t_uniq = p3i_output['reverse_tm'][j]
                amp_uniq = p3i_output['gold_amplicon'][j]
                poss_uniq = p3i_output['possible_off_targets'][j]
                amp_start_uniq = p3i_output['amplicon_start'][j]
                amp_end_uniq = p3i_output['amplicon_end'][j]
                amp_size_uniq = p3i_output['amplicon_size'][j]
                count_uniq = p3i_output['primer_count'][j]
                ideal_primer.append(ideal_uniq)
                alternate_primer.append(alt_uniq)
                forward_name.append(forward_n_uniq)
                forward_primer_list.append(forward_p_uniq)
                forward_tm_list.append(forward_t_uniq)
                reverse_name.append(reverse_n_uniq)
                reverse_primer_list.append(reverse_p_uniq)
                reverse_tm_list.append(reverse_t_uniq)
                amplicon_name.append(amp_uniq)
                off_target.append(poss_uniq)
                amplicon_start.append(amp_start_uniq)
                amplicon_end.append(amp_end_uniq)
                amplicon_size.append(amp_size_uniq)
                primer_count_high.append(count_uniq)
                continue

capture_best_matches()

def write_high_con():
    high_confidence_from_p3i = pd.DataFrame()
    high_confidence_from_p3i['amplicon'] = amplicon_name
    high_confidence_from_p3i['ideal (~195bp amplicon)'] = ideal_primer
    high_confidence_from_p3i['alternate (left or right)'] = alternate_primer
    high_confidence_from_p3i['forward_name'] = forward_name
    high_confidence_from_p3i['forward_primer'] = forward_primer_list
    high_confidence_from_p3i['forward_tm'] = forward_tm_list
    high_confidence_from_p3i['reverse_name'] = reverse_name
    high_confidence_from_p3i['reverse_primer'] = reverse_primer_list
    high_confidence_from_p3i['reverse_tm'] = reverse_tm_list
    high_confidence_from_p3i['off_targets'] = off_target
    high_confidence_from_p3i['amplicon_start'] = amplicon_start
    high_confidence_from_p3i['amplicon_end'] = amplicon_end
    high_confidence_from_p3i['amplicon_size'] = amplicon_size
    high_confidence_from_p3i['primer_count'] = primer_count_high
    high_confidence_from_p3i.to_csv(name+'_high_confidence_from_p3i.txt',sep='\t',index=False)
write_high_con()

print("Ideal finished proessing in --- %s seconds ---" % (time.time() - ideal_total_time))

def find_failed_primers_1():
    high_con = pd.read_csv(name+'_high_confidence_from_p3i.txt',delimiter='\t')
    high_con[['fam','chrom','pos','type']] = high_con['amplicon'].str.split('_',expand=True)
    high_con['variant_id'] = high_con['fam'] + '_' + high_con['chrom'] + '_' + high_con['pos']
    high_con.to_csv('high_confidence_with_variant_id',index=False)
    
    fam = []
    failed_chr = []
    failed_pos = []

    failed_1 = input_list[~input_list['variant_id'].isin(high_con['variant_id'])]

    failed_1.to_csv(name+'_failed_primers_right.csv',index=False)

find_failed_primers_1()

type_string = '_right'

def run_right_primers():

    primer_start_time = time.time()

    df = pd.read_csv(name+'_failed_primers_right.csv')
    add_600bp(df, reference)

    df.to_csv(name + '_601bp_region'+type_string+'.csv', index=False)

    write_primer3_file_right(df, name + '_primer3_input_right.txt',
                    60)

    with open(name + '_primer3_input_right.txt', 'rb', 0) as input_f, \
    open(name + '_human_output_right.txt', 'wb', 0) as output_f:
        subprocess.run(['primer3_core', '-format_output'], stdin=input_f,
                    stdout=output_f, check=True)

    with open(name + '_primer3_input_right.txt', 'rb', 0) as input_f, \
    open(name + '_machine_output_right.txt', 'wb', 0) as output_f:
        subprocess.run(['primer3_core'], stdin=input_f,
                    stdout=output_f, check=True)

    output = parse_p3_output(name + '_machine_output_right.txt', pairs=pairs)

    print("Right Primers Created in --- %s seconds ---" % (time.time() - primer_start_time))

    output.to_csv(name + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

run_right_primers()

row_split = []

split_primer_list()


#### Function below is for multiprocessing ###
#def ispcr_right(row):

#    print(row_split)
#    print(type_string)

#     print('Running isPcr')

#    start_time = time.time()

#    i = row

#    subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list_right.txt stdout -out=bed -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr_right.bed',shell=True)

#    subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list_right.txt stdout -out=fa -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr_right.fa',shell=True)

    #subprocess.run('isPcr GCF_000001405.40_GRCh38.p14_genomic.fna '+name+'_primer_list.txt stdout -out=bed -minPerfect=1 -minGood=6 -tileSize=11 -stepSize=5 >> '+name+'_isPcr_NC.bed',shell=True)
#ispcr_input = pd.read_csv(name + '_primer_list.txt',delimiter='\t',header=None)

#print(ispcr_input)

#ispcr(ispcr_input)

    #print("--- %s seconds ---" % (time.time() - start_time))
    
#pool = Pool(8)
#results = pool.map(ispcr_right, row_split, 8)

### End of function for Right-anchored multiprocessing ###

ispcr()

def process_right():

    right_total_time = time.time()

    ispcr_analysis()

    filter_ispcr_analysis()

    split_fasta()

    standard_present()

    standard_matches()

    percent_match()

    capture_best_matches()

    write_high_con()

    print("Right finished proessing in --- %s seconds ---" % (time.time() - right_total_time))

process_right()

def find_failed_primers_2():
    high_con = pd.read_csv(name+'_high_confidence_from_p3i.txt',delimiter='\t')
    high_con[['fam','chrom','pos','type']] = high_con['amplicon'].str.split('_',expand=True)
    high_con['variant_id'] = high_con['fam'] + '_' + high_con['chrom'] + '_' + high_con['pos']
    high_con.to_csv('high_confidence_with_variant_id',index=False)
    
    fam = []
    failed_chr = []
    failed_pos = []

    failed_2 = input_list[~input_list['variant_id'].isin(high_con['variant_id'])]

    failed_2.to_csv(name+'_failed_primers_left.csv',index=False)

find_failed_primers_2()

type_string = '_left'

def run_left_primers():

    primer_start_time = time.time()

    df = pd.read_csv(name+'_failed_primers_left.csv')
    add_600bp(df, reference)

    df.to_csv(name + '_601bp_region'+type_string+'.csv', index=False)

    write_primer3_file_left(df, name + '_primer3_input_left.txt',
                    60)

    with open(name + '_primer3_input_left.txt', 'rb', 0) as input_f, \
    open(name + '_human_output_left.txt', 'wb', 0) as output_f:
        subprocess.run(['primer3_core', '-format_output'], stdin=input_f,
                    stdout=output_f, check=True)

    with open(name + '_primer3_input_left.txt', 'rb', 0) as input_f, \
    open(name + '_machine_output_left.txt', 'wb', 0) as output_f:
        subprocess.run(['primer3_core'], stdin=input_f,
                    stdout=output_f, check=True)

    output = parse_p3_output(name + '_machine_output_left.txt', pairs=pairs)

    print("Left Primers Created in --- %s seconds ---" % (time.time() - primer_start_time))

    output.to_csv(name + '_primer_list'+type_string+'.txt', index=False, sep='\t',header=False)

run_left_primers()

row_split = []

split_primer_list()

### Function below is for Left-anchored multiprocessing ####
#def ispcr_left(row):
#def ispcr():

#    print(row_split)
#    print(type_string)


#    print('Running isPcr')

#    start_time = time.time()

    #for i in row_split:
    #for i in rows:

#    i = row

        #print('Running isPcr')

#    subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list_left.txt stdout -out=bed -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr_left.bed',shell=True)

#    subprocess.run('isPcr '+ref_path+' '+name+'_'+str(i)+'_primer_list_left.txt stdout -out=fa -minPerfect=1 -minGood=15 -tileSize=11 -stepSize=5 -maxSize=800 >> '+name+'_'+str(i)+'_isPcr_left.fa',shell=True)

    #subprocess.run('isPcr GCF_000001405.40_GRCh38.p14_genomic.fna '+name+'_primer_list.txt stdout -out=bed -minPerfect=1 -minGood=6 -tileSize=11 -stepSize=5 >> '+name+'_isPcr_NC.bed',shell=True)
#ispcr_input = pd.read_csv(name + '_primer_list.txt',delimiter='\t',header=None)

#print(ispcr_input)

#ispcr(ispcr_input)

#    print("--- %s seconds ---" % (time.time() - start_time))
#pool = Pool(8)
#results = pool.map(ispcr_left, row_split, 8)

#### End of function for Left-anchored multiprocessing ####

ispcr()

def process_left():

    left_total_time = time.time()

    ispcr_analysis()

    filter_ispcr_analysis()

    split_fasta()

    standard_present()

    standard_matches()

    percent_match()

    capture_best_matches()

    write_high_con()

    print("Left finished proessing in --- %s seconds ---" % (time.time() - left_total_time))

process_left()

high_confidence_from_p3i = pd.read_csv(name+'_high_confidence_from_p3i.txt', delimiter='\t')

high_confidence_from_p3i = high_confidence_from_p3i.drop_duplicates(subset=['amplicon']).reset_index()

primer_list_ideal = pd.read_csv(name+'_primer_list_ideal.txt',delimiter='\t', names = ['primer_name', 'forward','reverse'])

primer_list_right = pd.read_csv(name+'_primer_list_right.txt',delimiter='\t', names = ['primer_name', 'forward','reverse'])

primer_list_left = pd.read_csv(name+'_primer_list_left.txt',delimiter='\t', names = ['primer_name', 'forward','reverse'])


p3 = []

for i in range(len(input_list['variant_id'])):
    if ((input_list['variant_id'][i]+'_ideal') in high_confidence_from_p3i['amplicon'].values):
        p3_check = 'True'
        p3.append(p3_check)
        continue
    elif ((input_list['variant_id'][i]+'_right') in high_confidence_from_p3i['amplicon'].values):
        p3_check = 'True'
        p3.append(p3_check)
        continue
    elif ((input_list['variant_id'][i]+'_left') in high_confidence_from_p3i['amplicon'].values):
        p3_check = 'True'
        p3.append(p3_check)
        continue
    else:
        p3_check = 'False'
        p3.append(p3_check)
        continue

ispcr = []
primer = []
for i in range(len(input_list['variant_id'])):
    if ((input_list['variant_id'][i]+'_ideal') in high_confidence_from_p3i['amplicon'].values):
        p_name = input_list['variant_id'][i]+'_ideal'
        ispcr_check = 'True'
        ispcr.append(ispcr_check)
        primer.append(p_name)
        continue
    elif ((input_list['variant_id'][i]+'_right') in high_confidence_from_p3i['amplicon'].values):
        p_name = input_list['variant_id'][i]+'_right'
        ispcr_check = 'True'
        ispcr.append(ispcr_check)
        primer.append(p_name)
        continue
    elif ((input_list['variant_id'][i]+'_left') in high_confidence_from_p3i['amplicon'].values):
        p_name = input_list['variant_id'][i]+'_left'
        ispcr_check = 'True'
        ispcr.append(ispcr_check)
        primer.append(p_name)
        continue
    else:
        p_name = 'NA'
        ispcr_check = 'False'
        ispcr.append(ispcr_check)
        primer.append(p_name)
        continue

input_list['primer3'] = p3
input_list['isPcr'] = ispcr
input_list['primer_name'] = primer

ideal_primer_2 = []
alternate_primer_2 = []
forward_name_2 = []
forward_primer_list_2 = []
forward_tm_list_2 = []
reverse_name_2 = []
reverse_primer_list_2 = []
reverse_tm_list_2 = []
amplicon_name_2 = []
off_target_2 = []
amplicon_start_2 = []
amplicon_end_2 = []
amplicon_size_2 = []
primer_count_2 = []

for i in range(len(input_list['variant_id'])):
    if input_list['isPcr'][i] == 'False':
        ideal = 'NA_check_P3_primer_list'
        alt = 'NA'
        forward_n = 'NA'
        forward_p = 'NA'
        forward_t = 'NA'
        reverse_n = 'NA'
        reverse_p = 'NA'
        reverse_t = 'NA'
        amp = 'NA'
        poss_off = 'NA'
        amp_start = high_confidence_from_p3i['amplicon_start'][j]
        amp_end = high_confidence_from_p3i['amplicon_end'][j]
        amp_size  = 'NA'
        count = 'NA'
        ideal_primer_2.append(ideal)
        alternate_primer_2.append(alt)
        forward_name_2.append(forward_n)
        forward_primer_list_2.append(forward_p)
        forward_tm_list_2.append(forward_t)
        reverse_name_2.append(reverse_n)
        reverse_primer_list_2.append(reverse_p)
        reverse_tm_list_2.append(reverse_t)
        amplicon_name_2.append(amp)
        off_target_2.append(poss_off)
        amplicon_start_2.append(amp_start)
        amplicon_end_2.append(amp_end)
        amplicon_size_2.append(amp_size)
        primer_count_2.append(count)
    for j in range(len(high_confidence_from_p3i['amplicon'])):
        if input_list['variant_id'][i]+'_ideal' == high_confidence_from_p3i['amplicon'][j]:
                ideal = high_confidence_from_p3i['ideal (~195bp amplicon)'][j]
                alt = high_confidence_from_p3i['alternate (left or right)'][j]
                forward_n = high_confidence_from_p3i['forward_name'][j]
                forward_p = high_confidence_from_p3i['forward_primer'][j]
                forward_t = high_confidence_from_p3i['forward_tm'][j]
                reverse_n = high_confidence_from_p3i['reverse_name'][j]
                reverse_p = high_confidence_from_p3i['reverse_primer'][j]
                reverse_t = high_confidence_from_p3i['reverse_tm'][j]
                amp = high_confidence_from_p3i['amplicon'][j]
                poss_off = high_confidence_from_p3i['off_targets'][j]
                amp_start = high_confidence_from_p3i['amplicon_start'][j]
                amp_end = high_confidence_from_p3i['amplicon_end'][j]
                amp_size = high_confidence_from_p3i['amplicon_size'][j]
                count = high_confidence_from_p3i['primer_count'][j]
                ideal_primer_2.append(ideal)
                alternate_primer_2.append(alt)
                forward_name_2.append(forward_n)
                forward_primer_list_2.append(forward_p)
                forward_tm_list_2.append(forward_t)
                reverse_name_2.append(reverse_n)
                reverse_primer_list_2.append(reverse_p)
                reverse_tm_list_2.append(reverse_t)
                amplicon_name_2.append(amp)
                off_target_2.append(poss_off)
                amplicon_start_2.append(amp_start)
                amplicon_end_2.append(amp_end)
                amplicon_size_2.append(amp_size)
                primer_count_2.append(count)
                continue
        elif input_list['variant_id'][i]+'_ideal' not in high_confidence_from_p3i['amplicon'].values:
            if input_list['variant_id'][i]+'_right' == high_confidence_from_p3i['amplicon'][j]:
                ideal = high_confidence_from_p3i['ideal (~195bp amplicon)'][j]
                alt = high_confidence_from_p3i['alternate (left or right)'][j]
                forward_n = high_confidence_from_p3i['forward_name'][j]
                forward_p = high_confidence_from_p3i['forward_primer'][j]
                forward_t = high_confidence_from_p3i['forward_tm'][j]
                reverse_n = high_confidence_from_p3i['reverse_name'][j]
                reverse_p = high_confidence_from_p3i['reverse_primer'][j]
                reverse_t = high_confidence_from_p3i['reverse_tm'][j]
                amp = high_confidence_from_p3i['amplicon'][j]
                poss_off = high_confidence_from_p3i['off_targets'][j]
                amp_start = high_confidence_from_p3i['amplicon_start'][j]
                amp_end = high_confidence_from_p3i['amplicon_end'][j]
                amp_size = high_confidence_from_p3i['amplicon_size'][j]
                count = high_confidence_from_p3i['primer_count'][j]
                ideal_primer_2.append(ideal)
                alternate_primer_2.append(alt)
                forward_name_2.append(forward_n)
                forward_primer_list_2.append(forward_p)
                forward_tm_list_2.append(forward_t)
                reverse_name_2.append(reverse_n)
                reverse_primer_list_2.append(reverse_p)
                reverse_tm_list_2.append(reverse_t)
                amplicon_name_2.append(amp)
                off_target_2.append(poss_off)
                amplicon_start_2.append(amp_start)
                amplicon_end_2.append(amp_end)
                amplicon_size_2.append(amp_size)
                primer_count_2.append(count)
                continue
            elif (input_list['variant_id'][i]+'_right' not in high_confidence_from_p3i['amplicon'].values):
                if input_list['variant_id'][i]+'_left' == high_confidence_from_p3i['amplicon'][j]:
                    ideal = high_confidence_from_p3i['ideal (~195bp amplicon)'][j]
                    alt = high_confidence_from_p3i['alternate (left or right)'][j]
                    forward_n = high_confidence_from_p3i['forward_name'][j]
                    forward_p = high_confidence_from_p3i['forward_primer'][j]
                    forward_t = high_confidence_from_p3i['forward_tm'][j]
                    reverse_n = high_confidence_from_p3i['reverse_name'][j]
                    reverse_p = high_confidence_from_p3i['reverse_primer'][j]
                    reverse_t = high_confidence_from_p3i['reverse_tm'][j]
                    amp = high_confidence_from_p3i['amplicon'][j]
                    poss_off = high_confidence_from_p3i['off_targets'][j]
                    amp_start = high_confidence_from_p3i['amplicon_start'][j]
                    amp_end = high_confidence_from_p3i['amplicon_end'][j]
                    amp_size = high_confidence_from_p3i['amplicon_size'][j]
                    count = high_confidence_from_p3i['primer_count'][j]
                    ideal_primer_2.append(ideal)
                    alternate_primer_2.append(alt)
                    forward_name_2.append(forward_n)
                    forward_primer_list_2.append(forward_p)
                    forward_tm_list_2.append(forward_t)
                    reverse_name_2.append(reverse_n)
                    reverse_primer_list_2.append(reverse_p)
                    reverse_tm_list_2.append(reverse_t)
                    amplicon_name_2.append(amp)
                    off_target_2.append(poss_off)
                    amplicon_start_2.append(amp_start)
                    amplicon_end_2.append(amp_end)
                    amplicon_size_2.append(amp_size)
                    primer_count_2.append(count)
                    continue
        else:
                continue
print(len(forward_name_2))
print(len(p3))

input_list['ideal'] = ideal_primer_2
input_list['alternate'] = alternate_primer_2
input_list['forward_name'] = forward_name_2
input_list['forward_primer'] = forward_primer_list_2
input_list['forward_tm'] = forward_tm_list_2
input_list['reverse_name'] = reverse_name_2
input_list['reverse_primer'] = reverse_primer_list_2
input_list['reverse_tm'] = reverse_tm_list_2
input_list['amplicon_start'] = amplicon_start_2
input_list['amplicon_end'] = amplicon_end_2
input_list['amplicon_size'] = amplicon_size_2
input_list['primer_count'] = primer_count_2
input_list['possible_off_targets'] = off_target_2

input_list.to_csv(name+'_input_list_with_p3_ispcr_annotations.txt',sep='\t',index=False)

all_primers = pd.DataFrame()
all_primers['primers'] = all_primer_names
all_primers.to_csv('all_primers_list.csv',index=False)

def repseg_anno(df):

    print("Calculating Known Repetitive Region and Segdup Overlap Percentage")
    
    feat_bed = pd.DataFrame()
    feat_bed['chrom'] = df['CHROM']
    #feat_bed['chrom'] = feat_bed['chrom'].str.replace('chr','')
    feat_bed['pos1'] = df['POS'].astype(int)-1
    feat_bed['pos2'] = df['POS'].astype(int)-1
    
    feat_bed.to_csv('input_file_variant.bed',sep='\t',index=False,header=False)

    feat_bed_amp = pd.DataFrame()
    feat_bed_amp['chrom'] = df['CHROM']
    #feat_bed['chrom'] = feat_bed['chrom'].str.replace('chr','')
    feat_bed_amp['pos1'] = df['amplicon_start'].astype(int)-1
    feat_bed_amp['pos2'] = df['amplicon_end'].astype(int)

    feat_bed_amp.to_csv('input_file_amplicon.bed',sep='\t',index=False,header=False)

    subprocess.run('bedtools sort -i input_file_variant.bed > input_file_variant.sorted.bed', shell=True)

    subprocess.run('bedtools annotate -i input_file_variant.sorted.bed -files all_repeats.b38_new.bed segdup.hg38_new.bed > input_file_variant.repseg_anno.bed', shell=True)

    subprocess.run('bedtools sort -i input_file_amplicon.bed > input_file_amplicon.sorted.bed', shell=True)

    subprocess.run('bedtools annotate -i input_file_amplicon.sorted.bed -files all_repeats.b38_new.bed segdup.hg38_new.bed > input_file_amplicon.repseg_anno.bed', shell=True)

    primer_variant = pd.read_csv('input_file_variant.repseg_anno.bed',names=['chrom','pos1','pos2','all_repeat','segdup'],delimiter='\t')
    
    primer_variant['pos+'] = primer_variant['pos1']+1
    
    primer_variant['ID'] = primer_variant['chrom'].astype(str)+'_'+primer_variant['pos+'].astype(str)
    
    primer_variant_sort = primer_variant.sort_values(by=['ID'],ascending=True)

    primer_variant_sort.to_csv('primer_variant.txt',index=False,sep='\t')


    primer_amp = pd.read_csv('input_file_amplicon.repseg_anno.bed',names=['chrom','pos1','pos2','all_repeat','segdup'],delimiter='\t')
    
    primer_amp['pos+'] = primer_amp['pos1']+1
    
    primer_amp['ID'] = primer_amp['chrom'].astype(str)+'_'+primer_variant['pos+'].astype(str)
    
    primer_amp_sort = primer_amp.sort_values(by=['ID'],ascending=True)

    primer_amp_sort.to_csv('primer_amp.txt',index=False, sep='\t')

appended_input = pd.read_csv(name+'_input_list_with_p3_ispcr_annotations.txt',delimiter='\t')

repseg_anno(appended_input)

def merge_annotations():

    print("Annotating Primers with Known Repetitive Region and Segdup Overlap Percentage")
    
    df = pd.read_csv(name+'_input_list_with_p3_ispcr_annotations.txt',delimiter='\t')

    df['read_ID'] = df['CHROM'].astype(str)+'_'+df['POS'].astype(str)

    df2 = df.sort_values(by='read_ID', ascending=True)

    primer_variant_sort2 = pd.read_csv('primer_variant.txt',delimiter='\t')

    primer_amp_sort2 = pd.read_csv('primer_amp.txt',delimiter='\t')

    df2['variant_repeat'] = primer_variant_sort2['all_repeat'].values

    df2['variant_segdup'] = primer_variant_sort2['segdup'].values

    df2['variant_anno_ID'] = primer_variant_sort2['ID'].values

    df2['amplicon_repeat'] = primer_amp_sort2['all_repeat'].values

    df2['amplicon_segdup'] = primer_amp_sort2['segdup'].values

    df2['amplicon_anno_ID'] = primer_amp_sort2['ID'].values

    df2_sort = df2.sort_values(by=['CHROM','POS'])

    df2_sort.to_csv(name+'_FINAL_appended_input_file_with_all_anno.txt', sep='\t', index=False)

merge_annotations()



