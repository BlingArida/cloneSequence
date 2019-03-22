#!/usr/bin/env python3
'''
@Date: 2019-01-07 17:13:20
@LastEditors: cany
@LastEditTime: 2019-01-31 13:33:21

Extract CDR3 and sequence from change-o results
'''



import os,sys
from collections import defaultdict
import argparse


def get_index(alist, name):

    for idx, n in enumerate(alist):
        if n == name :
            return idx


def main():
    infile = args.get('infile') 
    outfile = args.get('outfile') or infile.rsplit('.', 1)[0]+'.cloneFraction.csv'
    Id = args.get('Id') 

    FractionList=[]
    dicCloneId = defaultdict(list)
    allLen = 0
    with open(infile,'r') as input:
        for line in input:
            lineinfo = line.strip().split('\t')
            if line.startswith('SEQUENCE_ID'):
                cloneId_idx = get_index(lineinfo,'CLONE')
                cdr3nc_idx = get_index(lineinfo,'JUNCTION')
                cdr3aa_idx = get_index(lineinfo,'JUNCTION_AA')
                cdr3aa_lenth = get_index(lineinfo,'JUNCTION_AA_LENGTH')
                headline = lineinfo             
                continue
            if lineinfo[cdr3aa_idx]!='':
                cloneInfo = lineinfo                
                #cloneInfo = [lineinfo[cdr3nc_idx],lineinfo[cdr3aa_idx],lineinfo[cdr3aa_lenth]]
                dicCloneId[lineinfo[cloneId_idx]].append(cloneInfo)                 


        
    for cloneId,line in dicCloneId.items():
        if Id:            
            if cloneId == Id:
                with open(infile.rsplit('.', 1)[0]+'cloneID_'+Id+'.tab','w') as idout:
                    idout.write('\t'.join(headline)+'\n')
                    for nline in line:
                        idout.write('\t'.join(nline)+'\n')            
        else:                
            cloneLen = len(line)        
            Fraction = [cloneId,cloneLen,line[0][cdr3nc_idx],line[0][cdr3aa_idx],line[0][cdr3aa_lenth]]####the first sequence in one clone group
            FractionList.append(Fraction)
            allLen += cloneLen

    if not Id:
        outlist = sorted(FractionList,key=lambda x:(x[1],x[0]),reverse=True)    
    
        with open(outfile,'w') as out:
            out.write('CloneId,CloneCounts,CloneFraction,CDR3_NC,CDR3_AA,CDR3_AA_lenth\n')
            for line in outlist:
                cloneFraction = format(float(line[1])/float(allLen),'.4f')
                line.insert(2,cloneFraction)
                stline = map(str,line)
                out.write(','.join(stline)+'\n')       
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser( 
        prog='clone_Fraction', 
        description='\t\tExtract CDR3 and Sequence', 
        epilog=( 
            'examples:\n  %(prog)s -i clone_inital.tab [output is clone_cluster file]\n' 
            '  %(prog)s -i clone_inital.tab -d cloneIdNumber [output is CDR3 and clone_sequence file]\n\n' 
            ), 
        formatter_class=argparse.RawTextHelpFormatter, 
    ) 
    parser.add_argument('-i', '--infile', help='The input tab file', required=True)
    parser.add_argument('-o', '--outfile', help='The output csv file',required = False) 
    parser.add_argument('-d', '--Id', help='The cloneId number',required=False) 

    args = vars(parser.parse_args()) 

    main() 
