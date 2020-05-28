#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 15:59:32 2018

@author: antoine
"""

from tqdm import *
import sys
import mmap

### FOR TQDM (Progress bar when reading large files)
def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def main():

    filec=sys.argv[1] # the file that contains cores 
    filep=sys.argv[2] # the file that contains peaks 
    
    COREs={} # dictonary that contains cores for each chromosome 
    newprep=[] #A list that contains all the peaks found in the cores
    
    #Store all COREs
 #Python checks if the specified key exists in the dict. If it does, then get() returns the value of that key. If the key does not exist, then get() returns the value specified in the second argument to get().
    with open(filec,"r") as core:
        for line in tqdm(core, total=get_num_lines(filec)):
            chrom,start,end = line.strip().split(" ")
            COREs[chrom]=COREs.get(chrom, [])
            COREs[chrom].append((start,end))
            
    #Store only peaks inside COREs        
    with open(filep,"r") as prep:
        for line in tqdm(prep, total=get_num_lines(filep)):
            chrom,start,end,pk,peak,score = line.strip().split("\t")
            if chrom in COREs:
                for (i,j) in COREs[chrom]:
                    if start>=i and end <=j:
                        newprep.append((chrom,start,end,pk,peak,score))
    
    #Write only peaks inside COREs                
    with open(filep.replace(".prepared",".fprepared"),"w") as file:
           for chrom,start,end,pk,peak,score in newprep:
                file.write(chrom+"\t"+start+"\t"+end+"\t"+pk+"\t"+peak+"\t"+score+"\n")
    
if __name__ == '__main__':
    main()
