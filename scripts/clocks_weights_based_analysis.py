import os, sys
from sklearn.preprocessing import StandardScaler
import math
import numpy as np
from scipy.stats import kstest
from scipy.stats import lognorm
from scipy.stats import mannwhitneyu
from scipy.stats import shapiro
from scipy import stats
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import json

try:
    if len(sys.argv) < 3:
        print("USAGE: python clocks_coefficient.txt path_to_intersected_data")
        exit(0)
except:
    pass
#
def process_coefficient_file(coefficient_file):
    fH = open (coefficient_file,"r")
    lines = fH.readlines()
    clocks = {}
    del lines[0]
    for line in lines:
        tokens = line.replace("\n","").split("\t")
        if tokens[0] != "Horvath-pan-mammalian":
            clocks[tokens[0]] = {"coefficients":[],"CpG":{}}
        else:
            clocks[tokens[0]+'-'+tokens[-1]] = {"coefficients":[],"CpG":{}}
    for line in lines:
        tokens = line.replace("\n","").split("\t")
        if tokens[0] != "Horvath-pan-mammalian":
            absolute_value = abs(float(tokens[2]))
            clocks[tokens[0]]["coefficients"].append(absolute_value)
            clocks[tokens[0]]["CpG"][tokens[1]] = {"coefficient":absolute_value ,"scaled_coefficient":"", "mutation":{"afr":"na","amr":"na","ami":"na","asj":"na","fin":"na","nfe":"na","eas":"na","sas":"na","mid":"na","sa":"na"}}
        else:
            absolute_value = abs(float(tokens[2]))
            clocks[tokens[0]+'-'+tokens[-1]]["coefficients"].append(absolute_value)
            clocks[tokens[0]+'-'+tokens[-1]]["CpG"][tokens[1]] = {"coefficient":absolute_value ,"scaled_coefficient":"", "mutation":{"afr":"na","amr":"na","ami":"na","asj":"na","fin":"na","nfe":"na","eas":"na","sas":"na","mid":"na","sa":"na"}}
    print("NOTICE: Done reading clocks file...")
    print("NOTICE: Found {} clocks...".format(len(clocks)))
    print("======================================")
    return clocks

#
def adding_mutation(path,clocks):
    exclude = ['Monika','Carola','Maria']
    files = os.listdir(path)
    for file in files:
        name = file.split('.')
        population = name[-4]
        fH = open (path+file,"r")
        lines = fH.readlines()
        if population != "sa":
            del lines[0]
        else:
            #continue
            print()
        for line in lines:
            tokens = line.replace("\n","").split(",")
            if tokens[1] not in exclude:
                    if tokens[1] == "Horvath-pan-mammalian":
                        for clock in clocks:
                            if "Horvath-pan-mammalian" in clock:
                                if tokens[2] in clocks[clock]["CpG"]:
                                    CpG_mutations = tokens[-1].split(";")
                                    if len(CpG_mutations) > 1:#multiple mutations...
                                        clocks[clock]["CpG"][tokens[2]]["mutation"][population] = "mutated"
                                    else:
                                        if CpG_mutations[0] != "": # only one mutation...
                                            clocks[clock]["CpG"][tokens[2]]["mutation"][population] = "mutated"
                                        else: # no mutation...
                                            clocks[clock]["CpG"][tokens[2]]["mutation"][population] = "not_mutated"
                    else:
                        CpG_mutations = tokens[-1].split(";")
                        if len(CpG_mutations) > 1:#multiple mutations...
                            clocks[tokens[1]]["CpG"][tokens[2]]["mutation"][population] = "mutated"
                        else:
                            if CpG_mutations[0] != "": # only one mutation...
                                clocks[tokens[1]]["CpG"][tokens[2]]["mutation"][population] = "mutated"
                            else: # no mutation...
                                clocks[tokens[1]]["CpG"][tokens[2]]["mutation"][population] = "not_mutated"
    print("NOTICE: Done adding the mutatations...")
    print("======================================")
    return clocks
#
def analyzing_clocks(clocks):
    fOut = open ("coefficients_stat.txt","w")
    fOut.write("clock\tmCpG_mean\tnmCpG_mean\tmCpG_median\tnmCpG_median\tnomalicy_test(shapiro)\tnon-parametric(mann-whitney u test)\n")

    data_scaled_all = {};data_all={}
    output = {}
    for clock in clocks:
        data_scaled_all[clock] = []
        data_all[clock]= {}
        output[clock] = {"scaled_data":{"normal_distribution_test":0,"enrichment_test":{"t-test":0,"non-parametric":0}} ,
                        "not_scaled_data":{"normal_distribution_test":0,"enrichment_test":{"t-test":0,"non-parametric":0}} 
                        }
        data = np.array(clocks[clock]["coefficients"])
        # Reshape the data to have two dimensions 
        #checked for the data order, and the ordered satyed unchanged...
        data = data.reshape(-1, 1)
        data_all[clock] = data
        #test for data distribution...
        result = kstest(data, 'norm')
        normalicy_stat,normalicy_p_value =  shapiro(data)
        # Initialize StandardScaler
        scaler = StandardScaler()
        # Fit scaler to your data to compute mean and standard deviation
        scaler.fit(data)
        # Transform your data using the fitted scaler
        data_scaled = scaler.transform(data)
        data_scaled_all[clock] = data_scaled
        #test for data distribution...
        result_scaled = kstest(data_scaled, 'norm')
        output[clock]["not_scaled_data"]["normal_distribution_test"] = result_scaled[1][0]
        # adding the scale data to clocks...
        for x in range(len(data)):
            for CpG in clocks[clock]["CpG"]:
                if data[x] == clocks[clock]["CpG"][CpG]["coefficient"]:
                    clocks[clock]["CpG"][CpG]["scaled_coefficient"] = data_scaled[x][0]
        print("NOTICE: Done scaling coefficients values...")
        #prepare data for stats...
        mCpG = []; nmCpG = []
        for CpG in clocks[clock]["CpG"]:
            # ckeck for mutations (in all the populations)...
            mutations_check = 0
            for population in clocks[clock]["CpG"][CpG]["mutation"]:
                if clocks[clock]["CpG"][CpG]["mutation"][population] == "mutated":
                    mutations_check +=1
            if mutations_check > 0:
                mCpG.append(clocks[clock]["CpG"][CpG]["coefficient"])
            else:
                nmCpG.append(clocks[clock]["CpG"][CpG]["coefficient"])
        #run stats...

        # non-parametric...
        non_parametric_p_value = mannwhitneyu(mCpG, nmCpG,alternative="two-sided").pvalue

        print("non-scaled nomalicy test:{}".format(normalicy_p_value))
        print("non-scaled non-parametric p-value:{}".format(non_parametric_p_value))
   
        #
        output[clock]["scaled_data"]["enrichment_test"]["non-parametric"] = non_parametric_p_value
        #output[clock]["not_scaled_data"]["enrichment_test"]["non-parametric"] = non_parametricp_value_scaled

        print("=================={}====================".format(clock))
  
        fOut.write("{}\t{}\t{}\t{}\t{}\t{}\t{}s\n".format(clock,np.mean(mCpG),np.mean(nmCpG),np.median(mCpG),np.median(nmCpG),normalicy_p_value,non_parametric_p_value))
        

    return 
#

if __name__ =="__main__":
    clocks = process_coefficient_file(sys.argv[1])
    clocks = adding_mutation(sys.argv[2],clocks)
    analyzing_clocks(clocks)
