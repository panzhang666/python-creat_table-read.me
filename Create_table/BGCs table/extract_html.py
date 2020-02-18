# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 09:30:05 2020

@author: mlk442
"""


import pandas as pd
import os
import re
import numpy as np
from functools import reduce
import copy
#from collection import defaultdict
import statistics
from bs4 import BeautifulSoup
from pathlib import Path
import os.path
import glob
import csv
from optparse import OptionParser
#extract digit for GenBank NODE_digit
def fist2digit(digit):
    if len(digit) >= 2:
        return (digit[0:2])
    else:
        return(digit+"_")
#add parameters
parser = OptionParser()
parser.add_option("-i", "--inputDir", dest="inputDir",
                  help="the folder for the outputs of antismash")
parser.add_option("-o", "--outputDir", dest="outdir",
                  help="the folder for the created BGC tables")
parser.add_option("-s", "--strain", dest="strain",
                  help="the antimash output of the specific strain")
parser.add_option("-g", "--gcfs_data", dest="path_gcf",
                  help="the antimash output of the specific strain")
(options, args) = parser.parse_args()
inputDir = options.inputDir
strain = options.strain
path_gcf = options.path_gcf

os.chdir(inputDir)
#os.chdir("C:\\Users\\mlk442\\Desktop\\manuscripts\\BGC_table\\")
#strain = "DA561A0323" 
#outdir = "C:\\Users\\mlk442\\Desktop\\manuscripts\\BGC_table\\"
outdir = options.outdir
csv_path =outdir +"/" +strain+"_BGC_table.csv"
with open (csv_path, "w") as f:
    csc_write = csv.writer(f)
    title ="BCGs in "+strain
    csv_head = ["BCGs in "+strain]
    csc_write.writerow(csv_head)

gbkfiles = []
for file in glob.glob(strain+"/*.gbk"):
    gbkfiles.append(file)
N_gbk = len(gbkfiles) -1
N_gbk      
#Read html file
path = strain +"/index.html"
htmlfile = open(path, "r", encoding = "utf-8")

#Creat gcfs(key) to gc(value) dictionary, which is specific to the strain
#Some gcs belong to more than one gcf
gcfsData = pd.read_table(path_gcf, sep = "\t")
dict_gcf = dict(zip(gcfsData.gcf.tolist(), gcfsData.bgcs.tolist()))
dict_gc = dict()
for key in dict_gcf.keys():
    gc_list = re.split(" ",dict_gcf[key])    
    filtered_list = list(filter(lambda x: re.match(strain, x) != None, gc_list))
    if len(filtered_list) >0 :
        for gc in filtered_list:
            if gc in dict_gc.keys():
               dict_gc[gc].append(key)
            else:
                dict_gc[gc] = [key]    
        
 #Process html file   
htmlhandle = htmlfile.read()
soup = BeautifulSoup(htmlhandle, "lxml")
count = 0
# Extract Contig Edge Warning information and make it a dictionary, Region(key)
dict_warn = dict()
warn = soup.find_all("div")
for i in range(len(warn)):
    if ("Region on contig edge." in warn[i].get_text()) and ("Download region GenBank file" in warn[i].get_text()):
        warn_input = re.split("\n",warn[i].get_text())
        #Four different kind format with the same information.
        #Here, I only use the one with "Region" in warn_input[1]
        if ("Region" in warn_input[1]):
            region_info = re.split(" |_|c",warn_input[1])
            name_region = "Region "+str(int(region_info[2]))+"."+region_info[7]
            dict_warn[name_region] = 0
            
dict_genbank = dict()
genbank = soup.find_all("div")
for i in range(len(genbank)):
    if ("Download region GenBank file" in genbank[i].get_text()):
        for link in genbank[i].find_all('a'):
            genbank_info = genbank[i].find_all('a')
            if len(genbank_info) >2:
                name_genbank = genbank_info[1].get("href")
                region_info2 = re.split("#|r|c",genbank_info[2].get("href"))
                name_region2 = "Region "+str(int(region_info2[2]))+"."+region_info2[3]
                dict_genbank[name_region2] = name_genbank 
                
df = pd.DataFrame(columns = ["Region","Type", "From", "To", "Contig Edge Warning","GenBank", "ID", "GCFs"])            
items = soup.find_all("tr")
faile=0
successed = 0
j=0
for i in range(len(items)):
    # Extract information from first 4 columns by "Region&nbsp" 
    if "Region&nbsp" in items[i].get_text(): 
        j += 1
        info = re.split("\n",items[i].get_text())
        Region = info[2].replace("&nbsp"," ")
        Type = info[5]
        From = info[7]
        To = info[8]
        #Find correaponding waring information in the warn dictionary
        if(Region in dict_warn.keys()):
            edgeWarn = "Region on contig edge"
        else:
            edgeWarn = ""
        #Creat GeneBank and ID information
        forGen = re.split(" |\.", Region)
        ID = strain +"_"+str(j)
        # Check whether Region exists or not, if not, add a row 
        if Region not in df.index.tolist():
            #Transform gdfs list to string
            GCFs = ", ".join(dict_gc[ID])
            GeneBank = strain +"/"+dict_genbank[Region]
            # Check if this Genbank file is in the antimash folder or not            
            if Path(GeneBank).is_file():
                df.loc[Region] = [Region,Type,From,To, edgeWarn, GeneBank, ID, GCFs]
                successed += 1
            else:
                faile += 0                 
#df.loc["BGCs"] = "BCGs in "+strain
#print (failed)
if (successed == N_gbk and faile ==0):
    print (strain +": completed successfully!")
else: 
    print (strain + ": "+str(N_gbk) + "gbk file in the folder, "+ str(success)+" gbk file in the table, " 
           + str(failed)+ " GenBank files can't be find in the directory!")    
df.to_csv(outdir +"/"+strain+"_BGC_table.csv", mode = "a", index = 0)
#get link


 
