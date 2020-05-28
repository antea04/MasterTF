

# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:28:56 2018
@author: antoine
"""

import sys
import os
import glob
import numpy as np
#import matplotlib.pyplot as plt
import pyBigWig
import time
import mmap
from tqdm import *
import re


### FOR TQDM (Progress bar when reading large files)
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
    interaction = sys.argv[5]
    core= open(sys.argv[6]) #/home/samar/NAS_Public/data/projects/Antoine_Networks/code/TFhubs_finder/data/pwm/core.txt
    hocomoco_threshold = open(sys.argv[7]) # /home/samar/NAS_Public/data/projects/Samira/SP1/HOCOMOCO_thresholds.txt to calculate the sarus score normalized 
    hoco_to_tf = open(sys.argv[8]) #/home/samar/NAS_Public/data/projects/Samira/MOUSE/HOCOMOCOv11_core_annotation_MOUSE_mono.tsv
    bigwigfile = sys.argv[9]
    genes_peaks_TFBS = {}
    peakchrom = {}
    TF_binding = {}
    listTF=[]
    peak_binding = {}
    gene_info = {}
    TF_occurence={}
    CORES={}
    realname={}
    thresholds ={}
    path=TFBSDir
    dictionary_overlapped_core = {}
    datapath=TFBSDir.replace("TFBS","data")
    #datapath = "/home/samar/NAS_Public/data/projects/Samira/mathematical_model/outputdir"
    for line in hoco_to_tf.readlines()[1:]:
        realname[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
    hoco_to_tf.close()


    for line in hocomoco_threshold.readlines()[1:]:
        name,p_value,_,_ = line.strip().split("\t")
        TF_name = realname[name]
        thresholds[TF_name] = p_value  

    hocomoco_threshold.close()

    filename=filep.split("/")[-1].replace(".promoters.bed", "")

    pb=float(0.003333333)

    delthresh=0.1
    
    if interaction =='True':
        delthresh=0.25
#def main():

#    filep=sys.argv[1]
#   prom = open(sys.argv[1])
 #   fant = open(sys.argv[2])
#  tads = open(sys.argv[3])
#    TFBSDir = sys.argv[4]
#    interaction = sys.argv[5]
#    
#   genes_peaks_TFBS = {}
#   peakchrom = {}
#   TF_binding = {}
#   listTF=[]
#   peak_binding = {}
 #   gene_info = {}
 #   path=TFBSDir
#  datapath=TFBSDir.replace("results/outputDir/TFBS","data")
#
#    filename=filep.split("/")[-1].replace(".promoters.bed", "")
    #pb=1/300
    #delthresh=0.1
    
    #if interaction =='True':
     #   delthresh=0.25
    #print(datapath)
    #print("delthresh=",delthresh)
    #print("reading gene peaks data... 1/8")

    #GETTING PROMOTER PEAKS INFO
    for line in prom.readlines():
        chrom,start,end,peak_name, gene_name = line.strip().split("\t")
        #gene_name = gene_name[:-1] #only for mouse, remove it for human 
        #print(line)
        #print(gene_name)
        for chunk in gene_name.split(";"):
            key, value = chunk.split("=")
            gene_info[key] = value
        gene = gene_info['gene_name']
        genes_peaks_TFBS[gene]=genes_peaks_TFBS.get(gene, {})  #gene -> peak_type,peaks
        genes_peaks_TFBS[gene]['prom']=genes_peaks_TFBS[gene].get('prom', set())
        genes_peaks_TFBS[gene]['prom'].add(peak_name)
        peakchrom[peak_name]=peakchrom.get(peak_name, {})  #peaks -> chrom,start,genes
        peakchrom[peak_name]['chrom']=chrom
        peakchrom[peak_name]['start']=start
        peakchrom[peak_name]['end']=end
        peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
        peakchrom[peak_name]['genes'].add(gene)

    #GETTING FANTOM5 ENHANCER PEAKS INFO
    for line in fant.readlines():
        chrom,start,end,peak_name, gene_name = line.strip().split("\t")
        gene_name = gene_name.split(";")
                
        if len(gene_name)>2: # it's two only for essay 
            gene = gene_name[2]        
        genes_peaks_TFBS[gene]=genes_peaks_TFBS.get(gene, {})
        genes_peaks_TFBS[gene]['fant']=genes_peaks_TFBS[gene].get('fant', set())
        genes_peaks_TFBS[gene]['fant'].add(peak_name)
        peakchrom[peak_name]=peakchrom.get(peak_name, {})
        peakchrom[peak_name]['chrom']=chrom
        peakchrom[peak_name]['start']=start
        peakchrom[peak_name]['end']=end
        peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
        peakchrom[peak_name]['genes'].add(gene)
    
    #GETTING TADS ENHANCER PEAKS INFO
    for line in tads.readlines():
        chrom,start,end,peak_name, gene_name = line.strip().split("\t")
        gene = gene_name.split("|")
        
        for i in gene:
            genes_peaks_TFBS[i]=genes_peaks_TFBS.get(i, {})
            genes_peaks_TFBS[i]['tads']=genes_peaks_TFBS[i].get('tads', set())
            genes_peaks_TFBS[i]['tads'].add(peak_name)
            peakchrom[peak_name]=peakchrom.get(peak_name, {})
            peakchrom[peak_name]['chrom']=chrom
            peakchrom[peak_name]['start']=start
            peakchrom[peak_name]['end']=end
            peakchrom[peak_name]['genes']=peakchrom[peak_name].get('genes', set())
            peakchrom[peak_name]['genes'].add(i)

    prom.close()
    fant.close()
    tads.close()
    for line in core.readlines():
        TF, s_core, e_core = line.strip().split("\t")
        CORES[TF]= (s_core, e_core)
    core.close() 
  
    print("reading TFBS data... 2/8")

    #bw = pyBigWig.open("/mnt/NAS_Public/data/projects/Samira/H1hES_validation/data/ATAC.H1hES.50Kcells.rep2.bwa.q20.rmdup_peaks_normalized.narrowPeak.wig.bw")
    bw = pyBigWig.open(bigwigfile)
    #bw = pyBigWig.open(os.path.abspath("/mnt/NAS_Public/data/projects/Antoine_Networks/data/" + filename + ".wig.bw"))  #BigWig file to get signal info for human 
    #bw= pyBigWig.open("/home/samar/NAS_Public/data/projects/Samira/SP1/data/01_3988Gustave_ATAC3-dox0_ATAC_mm_i11.merge200.pvalue005.bw") #for the mouse I am using 
    #Getting all motifs info from SARUS motif analysis
    #bw= pyBigWig.open("/home/samar/NAS_Public/data/projects/Samira/SP1/data/01_3988Gustave_ATAC3-dox0_ATAC_mm_i11.merge200.pvalue005.bw")
    for types in ["PROMOTERS","FANTOM5","TADS"]:
        print("from", types, "...")
        if types == "PROMOTERS":
            abvtypes = "prom"
        else:
            abvtypes = "enh"


        for tf in glob.glob(os.path.join(path+ "/" +types, "*.bed")):

            name = tf.split("/")[-1].replace(".bed","")

            
            if name in genes_peaks_TFBS:
                
                TFBS=open(tf)
                TF_binding[name]=TF_binding.get(name, {})
                for line in TFBS.readlines():

                    peak,pos1,pos2,_,score,sens = line.strip().split("\t")[:6]

                    #if score[0]!='-' and 1/(1+(1-pb)/(pb*2**float(score)))>delthresh:
                    #if score[0]!='-' and 1/(1+(1-pb)/(pb*2**float(score)))>delthresh: #1/(1+(1-pb)/(pb*2**float(score)))=probability of binding=affinity score
                        
                    peak_binding[peak]=peak_binding.get(peak, {})
                        
                        
                    peak_binding[peak][name]=peak_binding[peak].get(name, {})
                    peak_binding[peak][name][abvtypes]=peak_binding[peak][name].get(abvtypes, {})
                    peak_binding[peak][name][abvtypes][pos1+"-"+pos2]=peak_binding[peak][name][abvtypes].get(pos1+"-"+pos2, {})
                    peak_binding[peak][name][abvtypes][pos1+"-"+pos2]['affscore']=1/(1+(1-pb)/(pb*2**float(score)))
                    peak_binding[peak][name][abvtypes][pos1+"-"+pos2]['score_sarus']  = score
                    peak_binding[peak][name][abvtypes][pos1+"-"+pos2]['sens']=sens


                TFBS.close()
           
    r=0

    for peak in peak_binding:  #TF_occurence before filtering
        TF_occurence[peak]= TF_occurence.get(peak, {})
        for TF in peak_binding[peak]:
            TF_occurence[peak][TF]= TF_occurence[peak].get(TF, {})
            for types in peak_binding[peak][TF]:
                
                r= r + len(peak_binding[peak][TF][types])
            TF_occurence[peak][TF]['occu']= r
            r=0
    #listTF contains all actively transcribed TFs (more than one peak in the promoter before filtering)

    for i in TF_binding:
        if 'prom' in genes_peaks_TFBS[i]:
            if len(genes_peaks_TFBS[i]['prom'])>=1:   
                listTF.append(i)
    

    #middle of the peak, used to find position score
    with open(datapath + "/"+ filename, "r") as f:
    #with open("/home/samar/NAS_Public/data/projects/Samira/SP1/data/01_3988Gustave_ATAC3-dox0_ATAC_mm_i11.merge200.pvalue005_peaks.narrowPeak_threshold_kmeans.bed", "r" ) as f : #for mouse 
        for line in f.readlines():
            _,_,_,peak,_,_,_,_,_,mid = line.strip().split("\t")
            if peak in peak_binding:
                peak_binding[peak]['mid']=mid 

    dictionnary_to_compare = {}
    print(time.ctime())
    print("writing summary file")
    ### summary file read by networks.py
    
    #TF motifs are found in peaks, peaks are assigned to genes
    #We cross these 2 informations to get which TFs regulate which genes


    #Gsub_path=filep.replace("essai.bed",".matrix.summary.score_with_core_actively_transcribed_threshold1_essai.bed")
    Gsub_path=filep.replace("promoters.bed","resultat_final.bed")
    cpt=0
    tes=len(peak_binding)
    a="+"
    with open(Gsub_path,"w") as file:
        for peak in peak_binding:
            
            cpt+=1
            if cpt==int(0.25*tes):  #to have some progress check while writing
                print("25% ...")
            if cpt==int(0.5*tes):
                print("50% ...")
            if cpt==int(0.75*tes):
                print("75% ...")
                
            chrom=peakchrom[peak]['chrom']
            start=int(peakchrom[peak]['start'])
            mid=int(peak_binding[peak]['mid'])
            end = int(peakchrom[peak]['end'])
            
            
            for tf in peak_binding[peak]:
                if tf!='mid' and tf in listTF:

 #Only actively transcribed TFs are considerated
                    for types in peak_binding[peak][tf]:
 
                        for motif in peak_binding[peak][tf][types]:

                            pos1=int(motif.split("-")[0])
                            pos2=int(motif.split("-")[1])

                            genes=peakchrom[peak]['genes']

                            aff=peak_binding[peak][tf][types][motif]['affscore']
                            sarus_score = peak_binding[peak][tf][types][motif]['score_sarus']  

                            #affinity score
                            sarus_threshold = thresholds[tf]
                            strand = peak_binding[peak][tf][types][motif]['sens']
                            s_core = CORES[tf][0]
                            e_core = CORES[tf][1]
                            # I need the threshold too 
                            if str(strand) == a:
                                start_core = int(start)+ int(pos1) + int(s_core) -1 
                                end_core = int(start)+int(pos1) + int(e_core) -1
                            else: 
                                start_core = int(start)+int(pos2) - int(e_core) + 1
                                end_core = int(start)+int(pos2) - int(s_core) + 1
                            b=1
                            c=0

                            k=1
                            if bool(dictionnary_to_compare):
                                if peak in list(dictionnary_to_compare):
                                    if tf in list(dictionnary_to_compare[peak]):
                                        for t in list(dictionnary_to_compare[peak][tf]):

                                            for m in list(dictionnary_to_compare[peak][tf][t]):
                                                if b == 1: 


                                                    pos12 = int(m.split("-")[0])
                                                    pos22 = int(m.split("-")[1])
                                                    strand2 = dictionnary_to_compare[peak][tf][t][m]['sens']
                                                    if str(strand2) == a:
                                                        start_core2 = int(start)+ int(pos12) + int(s_core) -1 
                                                        end_core2 = int(start)+int(pos12) + int(e_core) -1
                                                    else: 
                                                        start_core2 = int(start)+int(pos22) - int(e_core) + 1
                                                        end_core2= int(start)+int(pos22) - int(s_core) + 1
 
                                #TF_grouped[((TF_grouped['end_core']<start)|(TF_grouped['start_core']> end))&(TF_grouped['strand']==x['strand'])] #filtering
                                                    if ((end_core2 < start_core) | (start_core2 > end_core)) : # the core of the new motif doesn't overlap with any other motif 
                                                        b=1
                                                        c=1
                                                    #elif: # the core of the new motif overlap with other motif
                                                    else:

                                                        if  (float(dictionnary_to_compare[peak][tf][t][m]['score_sarus']) < float(peak_binding[peak][tf][types][motif]['score_sarus'])): 

                                                            del dictionnary_to_compare[peak][tf][t][m]

                                                            b=1

                                                        else: 
                                                            b=0
                                                    if peak == "peak1" :
                                                        if tf == "BHLHE40" :
                                                            print("je suis peak", peak)
                                                            print("je suis tf", tf)
                                                            print("\n")
                                                            print("je suis peak_binding", peak_binding[peak][tf])
                                                            print("\n")
                                                            print("je suis dictionnay_to_compare", dictionnary_to_compare[peak][tf])
                                                            print("\n")
     #number_overlapped_TFs=number_overlapped_TFs+1
                            if b == 1:


                                dictionnary_to_compare[peak]=dictionnary_to_compare.get(peak, {})
                                
                                dictionnary_to_compare[peak][tf]=dictionnary_to_compare[peak].get(tf, {})
                                
                                dictionnary_to_compare[peak][tf][abvtypes]=dictionnary_to_compare[peak][tf].get(types, {}) #not inportant in filtering 
                                dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)]=dictionnary_to_compare[peak][tf][abvtypes].get(str(pos1)+"-"+str(pos2), {})
                                #dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)]= dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)].get(str(start_core)+"-"+str(end_core), {})
                                dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)]['affscore']=aff
                                dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)]['score_sarus']= sarus_score
                                dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)]['sens']=strand
                                dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)]= str(start_core)+"-"+str(end_core)

                    






                    for tipe in  dictionnary_to_compare[peak][tf]:

                        for mm in dictionnary_to_compare[peak][tf][tipe]:

                            number_overlapped_motifs = len(dictionnary_to_compare[peak][tf][tipe])

                            pos1=int(mm.split("-")[0])
                            pos2=int(mm.split("-")[1])

                            genes=peakchrom[peak]['genes']

                            aff=dictionnary_to_compare[peak][tf][tipe][mm]['affscore']
                            sarus_score = dictionnary_to_compare[peak][tf][tipe][mm]['score_sarus']  #affinity score
                            sarus_threshold = thresholds[tf]
                            strand = dictionnary_to_compare[peak][tf][tipe][mm]['sens']


                            for key in dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)].keys():
                                if (key != "sens") & (key != "affscore") & (key != "score_sarus"): 
                                    start_core = int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][key].split("-")[0])
                                    end_core = int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][key].split("-")[1])

                            #start_core = int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)].split("-")[0])
                            #end_core =int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)].split("-")[1])
                            number_overlapped_TFs=1
                            list_rank = [float(sarus_score)]
                            for TF in peak_binding[peak]: # the ranking part 
                                if TF != 'mid':
                                    if TF != tf:

                                        
                                        s_core1 = CORES[TF][0]
                                        e_core1 = CORES[TF][1]
                                        for t1 in peak_binding[peak][TF]:
                                            for m1 in peak_binding[peak][TF][t1]:
                                                pos13 = int(m1.split("-")[0])
                                                pos23 = int(m1.split("-")[1])
                                                strand3 = peak_binding[peak][TF][t1][m1]['sens']
                                                if str(strand3) == a:
                                                    start_core3 = int(start)+ int(pos13) + int(s_core1) -1 
                                                    end_core3 = int(start)+int(pos13) + int(e_core1) -1
                                                else: 
                                                    start_core3 = int(start)+int(pos23) - int(e_core1) + 1
                                                    end_core3= int(start)+int(pos23) - int(s_core1) + 1
                                            #sample[(sample['TF'] != name)&(sample['Chrom']==chrom)&(sample['end_core']>start)&(sample['start_core']<end)]
                                                if (end_core3 >= start_core) & (start_core3 <= end_core):
                                                    number_overlapped_TFs = number_overlapped_TFs + 1
                                                    sarus_score3 =  peak_binding[peak][TF][t1][m1]['score_sarus']
                                                    sarus_threshold3 = thresholds[TF]
                                                    list_rank.append(float(sarus_score3))


                                                    

                            seq= sorted(list_rank)
                            index = [seq.index(v) for v in list_rank]
                            rank = int(index[0])+1
                            print(chrom)
                            print(mid)
                            if chrom == "chrX" :
                                ch= "X"
                            elif chrom == "chrY":
                                ch = "Y"
                            else: 
                                ch = re.sub("[^0-9]", "", chrom)

                            if bw.stats(ch,start+pos1,start+pos2,exact=True)[0]==None:   #position score
                                pos=float(1/bw.values(ch,mid,mid+1)[0])      #this position score is higher if the motif is in the middle of the peak
                            else:                                                           #rather than on the edge of the peak
                            #for the narrowpeak I have from SH5Y5 cell type, the mid is already the position of the maximum height of the peak. Think to change 
                            #when dealing with other "true" narrowpeaks files to bw.values(ch, start+mid, start+mid+1)[0]

                                pos=float(bw.stats(ch,start+pos1,start+pos2,exact=True)[0]/bw.values(ch,mid,mid+1)[0])
                            score=pos*aff   #TFBS score
                            
                            minimum_intensity = min (bw.values(ch, start, start+1)[0], bw.values(ch, end, end+1)[0])
                            maximum_height = bw.values(ch, mid, mid+1)[0]
                            h = maximum_height - minimum_intensity
                            if peak== "peak1":
                                if TF == "BHLHE40":
                                    print ("je suis rank", rank)
                            file.write(chrom+ "\t"+tf+"\t"+str(start+pos1)+ "\t"+str(start+pos2)+"\t"+tf+"_"+','.join(genes)+"\t"+str((start+pos2) - (start+pos1))+"\t"+peak+"\t"+ str(end-start)+"\t"+str(pos)+"\t"+str(aff)+"\t"+str(sarus_score)+"\t"+str(sarus_threshold)+"\t"+str(maximum_height)+"\t"+str(h)+"\t"+str(bw.stats(ch,start+pos1,start+pos2,exact=True)[0])+"\t"+str(TF_occurence[peak][tf]['occu'])+"\t"+str(number_overlapped_motifs)+"\t"+str(number_overlapped_TFs)+"\t"+str(rank)+"\t"+str(strand)+"\t"+types+"\t"+ str(start_core)+"\t"+str(end_core)+"\n")


if __name__ == '__main__':
    main()