# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:02:17 2020

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

os.chdir("C:\\Users\\mlk442\\Desktop\\manuscripts\\correlation_table")
# Read tables
corTable = pd.read_table("cor_table.txt", sep = "\t")
gcfsData = pd.read_table("gcfs_data.txt", sep = "\t")
taxonomy = pd.read_table("taxonomy.txt", sep = "\t")
msTable = pd.read_table("ms_node_data.txt", sep = "\t", index_col = 0)

msNode_cor = corTable.groupby("ms_nodes")
msNodes = corTable.ms_nodes.unique().tolist()

dict_gcf = dict(zip(gcfsData.gcf.tolist(), gcfsData.bgcs.tolist()))
strainName = taxonomy.Genus +" "+ taxonomy.species
dict_strain = dict(zip(taxonomy.MaterialLabel.tolist(), strainName))

def pickRowname(rowname_list):
    rowname_list = reduce(lambda x,y: x+y, rowname_list)
#    rowname_list = rowname_list.split(" ")
    new_list = list()
    filtered_list = list(filter(lambda x: re.match("DA|DB", x) != None, rowname_list))
    for i in range(len(filtered_list)):
        new_list.append(filtered_list[i].split("_")[0])
    uniq_list= list(set(new_list))   
    return (uniq_list)

def extract_strain_from_spect(all_scan):
    new_scan = copy.deepcopy(all_scan)
    for i in range(len(all_scan)):
        new_scan[i] = re.split("-|_",all_scan[i])[1]
    return (new_scan)
    
def strain_spectra_mzs_intensity( all_scan, all_mz, all_intensity):
    # Creat dictionaries for spectra, m/z, intensity
    # Key is the strain, Values (mutiple) are the corresponding spectra, mz, and intensity  
    all_scan_strain = extract_strain_from_spect(all_scan)
    d_spectra = dict()
    d_mz = dict()
    d_intensity = dict()
    for i in range(len(all_scan_strain)):
        if all_scan_strain[i] in d_spectra.keys():
            d_spectra[all_scan_strain[i]].append(all_scan[i])
            d_mz[all_scan_strain[i]].append(all_mz[i])
            d_intensity[all_scan_strain[i]].append(all_intensity[i])
            
        else:
            d_spectra[all_scan_strain[i]]= [all_scan[i]]
            d_mz[all_scan_strain[i]]= [all_mz[i]]
            d_intensity[all_scan_strain[i]]= [all_intensity[i]]
            
    return(d_spectra, d_mz, d_intensity)
            

for ms in msNodes:
#    ms = "1023_DA977Z0659-SFM"
    msGroup = msNode_cor.get_group(ms)
#    ms = "1023_DA977Z0659-SFM"
    all_strain = msTable["strains"].loc[ms].split(" ")
    all_scan = msTable["all_scans"].loc[ms].split(" ")
#    all_scan_strain = extract_strain_from_spect(all_scan)
    all_mz = msTable["all_mzs"].loc[ms].split(" ")
    all_intensity = msTable["all_intensities"].loc[ms].split(" ") 
    dicts = strain_spectra_mzs_intensity(all_scan, all_mz, all_intensity)
    dict_spectra = dicts[0]
    dict_mz = dicts[1]
    dict_intensity = dicts[2]
 #   msGroup = msNode_cor.get_group(ms)
    gcfsName = msGroup.gcfs.tolist()
    # Rownames from GCFs cluster information
    rowname_list=list()
    for i in range(len(msGroup.gcfs)):
        addGC = dict_gcf[gcfsName[i]].split(" ")
        rowname_list.append(addGC)
        if (msGroup.gcf_knowns.tolist()[i] != "[]"): 
            gcfsName[i] = msGroup.gcfs.tolist()[i]+ " (" +msGroup.gcf_knowns.tolist()[i] +")" 
    rowname_list = pickRowname(rowname_list)
    if (not np.isnan(msGroup.ms_knowns.tolist()[1])):
        title = msGroup.ms_nodes+ "(" +msGroup.ms_knowns +") and its potential biosynthetic gene clusters:"
        title = title[0]
    else:
        title = ms +" (unknown) and its potential biosynthetic gene clusters:"
    # Add strain information from MS node information to rowname
    rowname_list=list(set(rowname_list+all_strain))
    colname = ["Species", "Spectra", "m/z (median)", "Ion Intensity (highest)"]
    colname= colname+gcfsName
    orig_gcfsName = msGroup.gcfs.tolist()
    #Cread dictionary for GCFs, key: GCFs correlated with specific msNodes, Value: BGC start with "DA or DB" 
    d_gcf = dict()
    for k in range(len(orig_gcfsName)):
        d_gcf[gcfsName[k]] = list(filter(lambda x: re.match("DA|DB", x) != None, 
                                    dict_gcf[orig_gcfsName[k]].split()))
    df = pd.DataFrame(columns = colname, index = rowname_list)
    for j in range(len(rowname_list)):
        df.iloc[j,0] = dict_strain[rowname_list[j]]
    #Fill in the blanks for GDFSs columns
        for h in range(len(gcfsName)):
            gc = list(filter(lambda x: re.match(rowname_list[j], x) != None, d_gcf[gcfsName[h]]))
            df.iloc[j,h+4] = " ".join(str(i) for i in gc)
        if (rowname_list[j] in dict_spectra.keys()):
            spect = dict_spectra[rowname_list[j]]
            df.iloc[j,1] = " ".join(str(i) for i in spect)           
            df.iloc[j,2] = statistics.median(list(map(eval,dict_mz[rowname_list[j]])))
            df.iloc[j,3] = max(list(map(eval, dict_intensity[rowname_list[j]])))
    corScore = msGroup.correlation_score.tolist()
    # sort by intensity
    df = df.sort_values(by = ["Ion Intensity (highest)"], ascending = False)
    df.loc["Metabologenomics correlation score (MS node - GCF)"]= ["","","",""]+corScore
    df.loc["MS node"] = title
    df.to_csv("cor_table.csv", mode = "a" )
    
    
    
            
            