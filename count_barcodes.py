#count_barcodes.py for mapping gRNA barcode

from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import glob
import datetime
import argparse

KEY = str.upper('gtttaaga')
reverse_read = False

def count_bcd(fastq_file, reference_file):
    print(datetime.datetime.now())
    print(f'fastq file: {fastq_file}')
    print(f'reference file: {reference_file}')
    
    #create files to write in
    ref = pd.read_excel(reference_file)
    ref_id = [i for i in ref['id']]
    ref_seq = [i.split('_')[1][-17:] for i in ref['id']]
    ref_dict = dict.fromkeys(ref_seq, 0)
    
    #open fastq file
    handle = open(fastq_file)
    readiter = SeqIO.parse(handle, "fastq")
    perfect_matches = 0
    non_perfect_matches = 0
    key_not_found = 0
    processed = 0
    
    print('\tcounting...\t\t', datetime.datetime.now())
    # process reads in fastq file
    for record in readiter:
        processed += 1
        if reverse_read:
            read = str.upper(str(record.seq.reverse_complement()))
        else:
            read = str.upper(str(record.seq))
        key_index = read.find(KEY) #if KEY not in read, -1 returned
        if key_index >= 0:
            start_index = key_index
            barcode = read[(start_index - 17):start_index]
            if barcode in ref_dict:
                ref_dict[barcode] += 1
                perfect_matches += 1
            else:
                non_perfect_matches += 1
        else:
            key_not_found += 1
            
        #if processed == 100000:
            #break
        
        
    print('\tcreating files...', '\t', datetime.datetime.now())
    # create csv file with guide id and respective counts in each row.
    df_count = pd.DataFrame(index = ref_id, columns = [fastq_file.split('.')[0]])
    for i in df_count.index:
        df_count.loc[i] = ref_dict[i.split('_')[1][-17:]]
    #df_count.index.name = 'sgRNA'
    df_count.to_csv('count_' + fastq_file.split('.')[0]+'.csv')

    # calculate statistical parameters
    percent_matched = round(perfect_matches/processed * 100, 2)
    percent_not_matched = round(non_perfect_matches/processed * 100, 2)
    percent_no_key = round(key_not_found/processed * 100, 2)
    guides_with_reads = np.count_nonzero(list(ref_dict.values()))
    guides_no_reads = len(ref_dict.values()) - guides_with_reads
    percent_no_reads = round(guides_no_reads/len(ref_dict.values()) * 100, 2)
        
    # write analysis statistics to csv file        
    df_st = pd.DataFrame(index = ['perfect matches','nonperfect match','key not found','processed reads',
                              'perfect matches (%)','nonperfect match (%)','key not found (%)','undetected guides (%)'])
    df_st[fastq_file.split('.')[0]] = [perfect_matches, non_perfect_matches, key_not_found, processed,
                                       percent_matched, percent_not_matched, percent_no_key, percent_no_reads]
    df_st.to_csv('statistics_' + fastq_file.split('.')[0]+'.csv')
    handle.close()
    print('\tfinised!\t\t', datetime.datetime.now(),'\n')
    return

path = r'' #dir of fastq and ref files
os.chdir(path)
ref = 'ref.xlsx'
fastq_files = glob.glob('*.fastq')

for f in fastq_files:
    count_bcd(f,ref)








