#!/usr/bin/env python3
'''
@LastEditors: cany
@LastEditTime: 2019-02-26 15:35:48
'''


import time
import numpy as np 
import pandas as pd
import os,sys,argparse
#from collections import Counter
import re



__author__ = 'yincan' 
__contact__ = 'yin_can@wuxiapptec.com' 



def get_index(alist, name):
    
    '''
    Get index of each column,using their title to locate
    '''

    for idx, n in enumerate(alist):
        if n == name :
            return idx


def ptm_cdr_wseq(seq,match_pattern):

    '''
    Cdrs
    AA such like "NG" "NS" "DG" on CDRs ;PTM means post translational modification.
    match_patten = r'N[GS]|DG'


    Wholeseq
    AA such like "N?S" "N?T" on the whole clone sequence.
    match_patten = r'N.[ST]'
    '''
    
    match = re.findall(match_pattern,seq)
    count = len(match)

    return count 


def cal_ptms(line):

    '''
    Calculate each parts PTM counts (first and second type).
    '''
    
    cdr1counts = ptm_cdr_wseq(line[0],r'N[GS]|DG')
    cdr2counts = ptm_cdr_wseq(line[1],r'N[GS]|DG')
    cdr3counts = ptm_cdr_wseq(line[2],r'N[GS]|DG')
    wseqcounts = ptm_cdr_wseq(line[3],r'N.[ST]')    
    
    return  cdr1counts,cdr2counts,cdr3counts,wseqcounts


def locC(df):

    '''
    There is supposed to be two CYS on a whole clone sequence. Locate the position of two CYS.
    '''

    ndf = (df == 'C').sum(axis = 0)
    loclist = ndf.sort_values(ascending = False).iloc[:2].index.tolist()
    
    return loclist
    

def indexC(row,loclist):

    '''
    Find each sequenc not the loc CYS location. Get ptm's count.
    '''
      
    clist = row[row == 'C'].index.tolist()    
    A = set(loclist) | set(clist)
    B = set(loclist) & set(clist)
    nlist = map(str,A - B)
    locate = ','.join(nlist)
    
    return len(A - B),locate


def main():
    infile = args.get('infile') 
    outfile = args.get('outfile') or infile.rsplit('.', 1)[0]+'.PTM.csv'

    
    t1 = time.time()
    clone_list = []

    
    with open(infile,'r') as input:
        for line in input:
            lineinfo = line.strip('').split('\t')
            if line.startswith('SEQUENCE_ID'):
                cdr1_idx = get_index(lineinfo,'CDR1_AA')
                cdr2_idx = get_index(lineinfo,'CDR2_AA')
                cdr3_idx = get_index(lineinfo,'JUNCTION_AA')
                cSeq_idx = get_index(lineinfo,'SEQUENCE_ALIGNMENT_AA')
                continue
            ### Discard sequening not cover all cdrsS
            if '' not in [lineinfo[cdr1_idx],lineinfo[cdr2_idx],lineinfo[cdr3_idx]]:
                cloneinfo = (lineinfo[cdr1_idx],lineinfo[cdr2_idx],lineinfo[cdr3_idx],lineinfo[cSeq_idx])
                clone_list.append(cloneinfo)

    ### Count the same number of sequences ,remove dup sequences to calculate less
    rmdup_clone = pd.value_counts(clone_list)
    #rmdup_clone = Counter(clone_list) 

    df_clone = pd.DataFrame({'cloneinfo':rmdup_clone.index,'Counts':rmdup_clone.values})  
    df_split_cloneinfo = pd.DataFrame((x for x in df_clone['cloneinfo']),index = df_clone.index,columns = ['Cdr1','Cdr2','Cdr3','Wseq'])
    df_merge_clone = pd.merge(df_split_cloneinfo,df_clone[['Counts']],right_index = True, left_index = True)

    
    ### Filter short clone seq,leave the longest sequence. 
    df_merge_clone['seqlen'] = df_merge_clone['Wseq'].map(lambda x : len(x))
    df_clone_filter = df_merge_clone[df_merge_clone['seqlen'] == df_merge_clone['seqlen'][0]]
   
    
    ### Locate C ,the max number of C counts index,largest C couny location ,supposed number = 2
    #df = df_clone_filter['Wseq'].str.split('',expand=True)
    df = pd.DataFrame(df_clone_filter['Wseq'].map(lambda x : list(x)).tolist())
    loclist = locC(df)
    print(loclist)
    df_clone_filter['Cptm'],df_clone_filter['Cptm_locate'] = zip(*df.apply(indexC,loclist = loclist,axis = 1))  
    
    df_clone_filter['Cdr1ptm'],df_clone_filter['Cdr2ptm'],df_clone_filter['Cdr3ptm'],df_clone_filter['Wseqptm'] \
    = zip(*df_clone_filter.apply(cal_ptms,axis=1))
    df_clone_filter['Allptms'] = df_clone_filter.Cptm + df_clone_filter.Cdr1ptm + df_clone_filter.Cdr2ptm \
    + df_clone_filter.Cdr3ptm + df_clone_filter.Wseqptm

    df_clone_filter.to_csv(outfile,index=0)
    t2 = time.time()

    print (t2-t1)     
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser( 
        prog='CalculatePTM', 
        description='\t\tFind PTMs and Sequence', 
        epilog=( 
            'examples:\n  %(prog)s -i clone_inital.tab [output is clone_cluster file]\n' 
            'contact: {} <{}>' 
        .format(__author__, __contact__)), 
        formatter_class=argparse.RawTextHelpFormatter, 
    ) 
    parser.add_argument('-i', '--infile', help='The input tab file', required=True)
    parser.add_argument('-o', '--outfile', help='The output csv file',required = False) 

    args = vars(parser.parse_args()) 

    main() 
