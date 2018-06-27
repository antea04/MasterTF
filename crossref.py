#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:28:56 2018

@author: antoine
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
import time
import mmap
from tqdm import *



def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def main():

    filep=sys.argv[1]
    prom = open(sys.argv[1])
    fant = open(sys.argv[2])
    tads = open(sys.argv[3])
    TFBSDir = sys.argv[4]
    #filep="../../results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.promoters.bed"
    #prom = open("../../results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.promoters.bed")
    #fant = open("../../results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.FANTOM5_enhanc.bed")
    #tads = open("../../results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.TADS_enhanc.bed")
    #TFBSDir = "../../results/outputDir/TFBS/"
    genes_peaks_TFBS = {}
    peakchrom = {}
    TF_binding = {}
    peak_binding = {}
    gene_info = {}
    path=TFBSDir
    filepath=TFBSDir.replace("results/outputDir/TFBS","data")
    filename=filep.split("/")[-1].replace(".promoters.bed", "")
    pb=1/300
    delthresh=0.2
    print("delthresh=",delthresh)
    print("reading gene peaks data... 1/8")
    
    for line in prom.readlines():
        chrom,start,_,peak_name, gene_name = line.strip().split("\t")
        for chunk in gene_name.split(";"):
            key, value = chunk.split("=")
            gene_info[key] = value
        gene = gene_info['gene_name']
        genes_peaks_TFBS[gene]=genes_peaks_TFBS.get(gene, {})
        genes_peaks_TFBS[gene]['prom']=genes_peaks_TFBS[gene].get('prom', set())
        genes_peaks_TFBS[gene]['prom'].add(peak_name)
        peakchrom[peak_name]=peakchrom.get(peak_name, {})
        peakchrom[peak_name]['chrom']=chrom
        peakchrom[peak_name]['start']=start
        peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
        peakchrom[peak_name]['genes'].add(gene)
    
    for line in fant.readlines():
        chrom,start,_,peak_name, gene_name = line.strip().split("\t")
        gene_name = gene_name.split(";")
        if len(gene_name)>2:
            gene = gene_name[2]        
        genes_peaks_TFBS[gene]=genes_peaks_TFBS.get(gene, {})
        genes_peaks_TFBS[gene]['enh']=genes_peaks_TFBS[gene].get('enh', set())
        genes_peaks_TFBS[gene]['enh'].add(peak_name)
        peakchrom[peak_name]=peakchrom.get(peak_name, {})
        peakchrom[peak_name]['chrom']=chrom
        peakchrom[peak_name]['start']=start
        peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
        peakchrom[peak_name]['genes'].add(gene)
    
    for line in tads.readlines():
        chrom,start,_,peak_name, gene_name = line.strip().split("\t")
        gene = gene_name.split("|")
        for i in gene:   
            genes_peaks_TFBS[i]=genes_peaks_TFBS.get(i, {})
            genes_peaks_TFBS[i]['enh']=genes_peaks_TFBS[i].get('enh', set())
            genes_peaks_TFBS[i]['enh'].add(peak_name)
            peakchrom[peak_name]=peakchrom.get(peak_name, {})
            peakchrom[peak_name]['chrom']=chrom
            peakchrom[peak_name]['start']=start
            peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
            peakchrom[peak_name]['genes'].add(i)
    
    
    prom.close()
    fant.close()
    tads.close()           
    
    print("reading TFBS data... 2/8")
    
    bw = pyBigWig.open(os.path.abspath(filepath+filename+".wig.bw"))
    
    for types in ["promoters","FANTOM5","TADS"]:
        print("from", types, "...")
        if types == "promoters":
            abvtypes = "prom"
        else:
            abvtypes = "enh"
        for tf in glob.glob(os.path.join(path+types, "*.bed")):
            if tf.split("/")[-1].replace(".bed","")=="NANOG2":
                name = "NANOG"
            else:
                name = tf.split("/")[-1].replace(".bed","")
            if name in genes_peaks_TFBS:
                TFBS=open(tf)
                TF_binding[name]=TF_binding.get(name, {})
                TF_binding[name][abvtypes]=TF_binding[name].get(abvtypes, {})
                for line in TFBS.readlines():
                    peak,pos1,pos2,_,score,sens = line.strip().split("\t")[:6]
                    if score[0]!='-' and 1/(1+(1-pb)/(pb*2**float(score)))>delthresh:
                        peak_binding[peak]=peak_binding.get(peak, {})
                        peak_binding[peak][name]=peak_binding[peak].get(name, {})
                        peak_binding[peak][name][abvtypes]=peak_binding[peak][name].get(abvtypes, {})
                        peak_binding[peak][name][abvtypes][pos1+"-"+pos2]=peak_binding[peak][name][abvtypes].get(pos1+"-"+pos2, {})
                        peak_binding[peak][name][abvtypes][pos1+"-"+pos2]['affscore']=1/(1+(1-pb)/(pb*2**float(score)))
                        peak_binding[peak][name][abvtypes][pos1+"-"+pos2]['sens']=sens
                TFBS.close()
    
    with open(filepath+filename, "r") as f:  
        for line in f.readlines():
            _,_,_,peak,_,_,_,_,_,mid = line.strip().split("\t")
            if peak in peak_binding:
                peak_binding[peak]['mid']=mid 
    
    print(time.ctime())
    print("writing cross file")
    Gsub_path=filep.replace(".promoters.bed",".cross")
    cpt=0
    tes=len(peak_binding)
#    listprob=[]
#    listposi=[]
    with open(Gsub_path+".bed","w") as file:
        for peak in peak_binding:
            cpt+=1
            if cpt==int(0.25*tes):
                print("25% ...")
            if cpt==int(0.5*tes):
                print("50% ...")
            if cpt==int(0.75*tes):
                print("75% ...")
            chrom=peakchrom[peak]['chrom']
            start=int(peakchrom[peak]['start'])
            mid=int(peak_binding[peak]['mid'])
            for tf in peak_binding[peak]:
                if tf!='mid':
                    for types in peak_binding[peak][tf]:
                        for motif in peak_binding[peak][tf][types]:
                            pos1=int(motif.split("-")[0])
                            pos2=int(motif.split("-")[1])
                            genes=peakchrom[peak]['genes']
                            aff=peak_binding[peak][tf][types][motif]['affscore']
#                            listprob.append(aff)
                            pos=float(bw.stats(chrom,start+pos1,start+pos2,exact=True)[0]/bw.values(chrom,start+mid,start+mid+1)[0])
#                            listposi.append(pos)
                            score=pos*aff
                            sens=peak_binding[peak][tf][types][motif]['sens']
                            if score >= 0.1:
                                file.write(chrom+"\t"+str(start+pos1)+"\t"+str(start+pos2)+"\t"+tf+"_"+','.join(genes)+"\t"+str(score)+"\t"+sens+"\t"+peak+"\t"+types+"\n")
   
#    plt.figure()
#    plt.hist(listprob, bins=20, label="Binding Probability")
#    plt.show(block=False)
#    plt.figure()
#    plt.hist(listposi, bins=20, label="Position Score")
#    plt.show()

    
if __name__ == '__main__':
    main()