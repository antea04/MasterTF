#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import glob

def main():
    
    thresholds_input = open(sys.argv[1])
    regions_type = sys.argv[4]
    input_fasta = open(sys.argv[3])

    pwmDir = sys.argv[2]
    thresholds_input = open(sys.argv[1])
    hoco_to_tf = open(sys.argv[5])
    thresholds = {}
    realname = {}
    #Gsub_path2 = input_fasta.name + "corrected"
    #with open(Gsub_path2,"w") as file:
        #for line in input_fasta.readlines():
            #if line.strip().split(":"):
                #peak = line.strip().split(":")[0]
                #file.write(peak+ "\n")
            #else:
                #if not (":" in line):
                    #file.write(line+"\n")
    for line in thresholds_input.readlines()[1:]:
        name,_,p_value,_ = line.strip().split("\t")
        thresholds[name] = p_value  
    thresholds_input.close()
    path = pwmDir
    
    for line in hoco_to_tf.readlines()[1:]:
        realname[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
    hoco_to_tf.close()
    
    scriptPath=os.path.dirname(os.path.realpath(__file__))

    outPath=os.path.dirname(os.path.realpath(input_fasta.name))
    outPath = outPath.replace('FASTA', 'TFBS')
    outPath = outPath+"/"+regions_type+"/"
    print(outPath)
    if not os.path.exists(outPath):
        os.mkdir(outPath)

    for tf in glob.glob(os.path.join(path, "*.pwm")):
        
        name = tf.split("/")[-1].replace(".pwm","")
        
        threshold = thresholds[name]
        
        command='java -cp '+ scriptPath+'sarus-01Mar2018.jar ru.autosome.SARUS '+ input_fasta.name +' '+ tf+' '+ threshold+ ' --output-bed skipn'
        subprocess.call(["java", "-cp", scriptPath+"/sarus-01Mar2018.jar", "ru.autosome.SARUS", input_fasta.name, tf, threshold, "--output-bed", "skipn"], stdout=open(outPath+realname[name]+".bed", "w"))
        #print(str(command))
        #subprocess.call(["java", "-cp", scriptPath+"/sarus-01Mar2018.jar", "ru.autosome.SARUS", input_fasta.name, tf, threshold, "--output-bed", "skipn"], stdout=open(outPath+realname[name]+".bed", "w"))

if __name__ == '__main__':
    main()

