# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:28:56 2018; modified on June 4, 2020
@author: Antoine and Valentina
"""

import sys
import os
import glob
import numpy as np
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

    fileWithPeaksBed6 = sys.argv[1] #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/testOutput/ATAC.GM12878.50Kcells.rep1_peaks.narrowPeak.bed6.bed
    TFBSDir = sys.argv[2] #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/testOutput/TFBS
    core= open(sys.argv[3]) #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/annotationFiles/core.txt
    hocomoco_threshold = open(sys.argv[4]) # /home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/annotationFiles/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt to calculate the sarus score normalized 
    hoco_to_tf = open(sys.argv[5]) #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/annotationFiles/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv
    bigwigfile = sys.argv[6] #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/testData/ATAC.GM12878.50Kcells.rep1.wig.bw
    expressedGenes = sys.argv[7] #/home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject/testOutput/peaksInExpressedGenes.bed
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
    datapath=TFBSDir.replace("TFBS","")
    for line in hoco_to_tf.readlines()[1:]:
        realname[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
    hoco_to_tf.close()


    for line in hocomoco_threshold.readlines()[1:]:
        name,p_value,_,_ = line.strip().split("\t")
        TF_name = realname[name]
        thresholds[TF_name] = p_value  

    hocomoco_threshold.close()

    pb=float(0.003333333)
   

    #GETTING PEAKS INFO
    
    with open(fileWithPeaksBed6, "r") as f:
        for line in f.readlines():    
            chrom,start,end,peak_name = line.strip().split("\t")[0:4]
            peakchrom[peak_name]=peakchrom.get(peak_name, {})  #peaks -> chrom,start,genes
            peakchrom[peak_name]['chrom']=chrom
            peakchrom[peak_name]['start']=start
            peakchrom[peak_name]['end']=end
        f.close()

    #GETTING COORDINATES OF MOTIF CORES
    for line in core.readlines():
        TF, s_core, e_core = line.strip().split("\t")
        CORES[TF]= (s_core, e_core)
    core.close() 
  
    
  
    bw = pyBigWig.open(bigwigfile)
    
    
    
    #getting expressed genes from expressedGenes
    with open(expressedGenes, "r") as f:
        for line in f.readlines():
            _,_,_,_,geneName = line.strip().split("\t")
            if (geneName not in listTF):
                listTF.append(geneName)        
        f.close()
    
    
 
    
    
    #Getting all motifs info from SARUS motif analysis

    for tf in glob.glob(os.path.join(path+ "/" +"regions", "*.bed")):

        name = tf.split("/")[-1].replace(".bed","")  

        if name in listTF:

            TFBS=open(tf)
            TF_binding[name]=TF_binding.get(name, {})
            for line in TFBS.readlines():
                peak,pos1,pos2,_,score,sens = line.strip().split("\t")[:6]                      
                peak_binding[peak]=peak_binding.get(peak, {})
                peak_binding[peak][name]=peak_binding[peak].get(name, {})
                peak_binding[peak][name][pos1+"-"+pos2]=peak_binding[peak][name].get(pos1+"-"+pos2, {})
                peak_binding[peak][name][pos1+"-"+pos2]['affscore']=1/(1+(1-pb)/(pb*2**float(score)))
                peak_binding[peak][name][pos1+"-"+pos2]['score_sarus']  = score
                peak_binding[peak][name][pos1+"-"+pos2]['sens']=sens
            TFBS.close()
        else:
            print ("will ignore TF: {}".format(name))
           

    for peak in peak_binding:  #TF_occurence before filtering
        TF_occurence[peak]= TF_occurence.get(peak, {})
        for TF in peak_binding[peak]:
            TF_occurence[peak][TF]= TF_occurence[peak].get(TF, {})                
            TF_occurence[peak][TF] = len(peak_binding[peak][TF])
     

    #middle of the peak, used to find position score
    with open(fileWithPeaksBed6, "r") as f:
        for line in f.readlines():
            _,_,_,peak,_,mid = line.strip().split("\t")
            if peak in peak_binding:
                peak_binding[peak]['mid']=mid 
        f.close()

    print(time.ctime())
    
    #TF motifs are found in peaks, peaks are assigned to genes
    #We cross these 2 informations to get which TFs regulate which genes


    Gsub_path=fileWithPeaksBed6.replace("bed6.bed","final_result.bed")
    cpt=0
    tes=len(peak_binding)
    a="+"
    with open(Gsub_path,"w") as file:
        
        #print the header:
            
        file.write("chrom\tstartMotif\tendMotif\tTF\tscore\tstrand\tmotifLength\tpeak\t"+ "peakLenth\tpositionalScore\tprobabilityAffinityScore\tSarusAffinityScore\tSarusThreshold\tpeakMaximumHeight\t"+ "peakHeightDelta\tmotifSignalHeight\tnumberOfTFHitsPerPeak\tnumberOfNonOverlappingMotifs\tnumberOverlappingTFs\trankAmongOverlappingMotifs\tmotifCoreStart\tmotifCoreEnd\n")

        
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
                
                
                dictionnary_to_compare = {}
    
                
                for tf in peak_binding[peak]:
                    if tf!='mid' and tf in listTF:
    
                        #Only actively transcribed TFs are considerated
     
                        for motif in peak_binding[peak][tf]:
    
                            pos1=int(motif.split("-")[0])
                            pos2=int(motif.split("-")[1])
    
                            aff=peak_binding[peak][tf][motif]['affscore']
                            sarus_score = peak_binding[peak][tf][motif]['score_sarus']  
    
                            #affinity score
                            sarus_threshold = thresholds[tf]
                            strand = peak_binding[peak][tf][motif]['sens']
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
                                if tf in list(dictionnary_to_compare):
                                    for m in list(dictionnary_to_compare[tf]):
                                            
                                        if b == 1: 
    
                                            pos12 = int(m.split("-")[0])
                                            pos22 = int(m.split("-")[1])
                                            strand2 = dictionnary_to_compare[tf][m]['sens']
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
    
                                                if  (float(dictionnary_to_compare[tf][m]['score_sarus']) < float(peak_binding[peak][tf][motif]['score_sarus'])): 
    
                                                    del dictionnary_to_compare[tf][m]
    
                                                    b=1
    
                                                else: 
                                                    b=0
                                   
                                    #number_overlapped_TFs=number_overlapped_TFs+1
                            if b == 1:
    
                                dictionnary_to_compare[tf]=dictionnary_to_compare.get(tf, {})
                                    
                                dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)]=dictionnary_to_compare[tf].get(str(pos1)+"-"+str(pos2), {})
                                    #dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)]= dictionnary_to_compare[peak][tf][abvtypes][str(pos1)+"-"+str(pos2)].get(str(start_core)+"-"+str(end_core), {})
                                dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)]['affscore']=aff
                                dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)]['score_sarus']=sarus_score
                                dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)]['sens']=strand
                                dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)]= str(start_core)+"-"+str(end_core)
    
    
    
                        for mm in dictionnary_to_compare[tf]:
    
                            number_overlapped_motifs = len(dictionnary_to_compare[tf])
    
                            pos1=int(mm.split("-")[0])
                            pos2=int(mm.split("-")[1])
    
                            aff=dictionnary_to_compare[tf][mm]['affscore']
                            sarus_score = dictionnary_to_compare[tf][mm]['score_sarus']  #affinity score
                            sarus_threshold = thresholds[tf]
                            strand = dictionnary_to_compare[tf][mm]['sens']
    
    
                            for key in dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)].keys():
                                if (key != "sens") & (key != "affscore") & (key != "score_sarus"): 
                                    start_core = int(dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)][key].split("-")[0])
                                    end_core = int(dictionnary_to_compare[tf][str(pos1)+"-"+str(pos2)][key].split("-")[1])
    
                            #start_core = int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)].split("-")[0])
                            #end_core =int(dictionnary_to_compare[peak][tf][tipe][str(pos1)+"-"+str(pos2)][str(start_core) +"-"+str(end_core)].split("-")[1])
                            number_overlapped_TFs=1
                            list_rank = [float(sarus_score)]
                            for TF in peak_binding[peak]: # the ranking part 
                                if TF != 'mid':
                                    if TF != tf:
                                            
                                        s_core1 = CORES[TF][0]
                                        e_core1 = CORES[TF][1]
                                        for m1 in peak_binding[peak][TF]:
                                            pos13 = int(m1.split("-")[0])
                                            pos23 = int(m1.split("-")[1])
                                            strand3 = peak_binding[peak][TF][m1]['sens']
                                            if str(strand3) == a:
                                                start_core3 = int(start)+ int(pos13) + int(s_core1) -1 
                                                end_core3 = int(start)+int(pos13) + int(e_core1) -1
                                            else: 
                                                start_core3 = int(start)+int(pos23) - int(e_core1) + 1
                                                end_core3= int(start)+int(pos23) - int(s_core1) + 1
                                            #sample[(sample['TF'] != name)&(sample['Chrom']==chrom)&(sample['end_core']>start)&(sample['start_core']<end)]
                                            if (end_core3 >= start_core) & (start_core3 <= end_core):
                                                number_overlapped_TFs = number_overlapped_TFs + 1
                                                sarus_score3 =  peak_binding[peak][TF][m1]['score_sarus']
                                                sarus_threshold3 = thresholds[TF]
                                                list_rank.append(float(sarus_score3))
                                                    
                            seq= sorted(list_rank,  reverse=True)
                            index = [seq.index(v) for v in list_rank]
                            rank = int(index[0])+1
    
                            if chrom == "chrX" :
                                ch= "X"
                            elif chrom == "chrY":
                                ch = "Y"
                            else: 
                                ch = re.sub("[^0-9]", "", chrom)
                            if ch not in bw.chroms():
                                ch=chrom
    
                            if bw.stats(ch,start+pos1,start+pos2,exact=True)[0]==None:   #position score
                                pos=float(1/bw.values(ch,start+mid,start+mid+1)[0])      #this position score is higher if the motif is in the middle of the peak
                            else:                                                           #rather than on the edge of the peak
                                #for the narrowpeak I have from SH5Y5 cell type, the mid is already the position of the maximum height of the peak. Think to change 
                                #when dealing with other "true" narrowpeaks files to bw.values(ch, start+mid, start+mid+1)[0] - VALENTINA: DID THIS CHANGE!!
    
                                pos=float(bw.stats(ch,start+pos1,start+pos2,exact=True)[0]/bw.values(ch,start+mid,start+mid+1)[0])
                            score=pos*aff   #TFBS score
                                
                            minimum_intensity = min (bw.values(ch, start, start+1)[0], bw.values(ch, end, end+1)[0])
                            maximum_height = bw.values(ch, start+ mid,start + mid+1)[0] #VALENTINA: changed Samira's code here
                            h = maximum_height - minimum_intensity
     
                            file.write(chrom +"\t"+str(start+pos1)+ "\t"+str(start+pos2)+ "\t"+tf + "\t" + str(score)+"\t"+str(strand) +"\t"+str((start+pos2) - (start+pos1))+"\t"+peak+"\t"+ str(end-start)+"\t"+str(pos)+"\t"+str(aff)+"\t"+str(sarus_score)+"\t"+str(sarus_threshold)+"\t"+str(maximum_height)+"\t"+str(h)+"\t"+str(bw.stats(ch,start+pos1,start+pos2,exact=True)[0])+"\t"+str(TF_occurence[peak][tf])+"\t"+str(number_overlapped_motifs)+"\t"+str(number_overlapped_TFs)+"\t"+str(rank)+"\t"+ str(start_core)+"\t"+str(end_core)+"\n")


if __name__ == '__main__':
    main()
