
import os, sys
import argparse
import sklearn
import numpy as np
import pandas as pd
from biolearn.data_library import DataLibrary
from biolearn.model_gallery import ModelGallery
import seaborn as sn
import urllib.request
from urllib.request import urlopen
import ssl
import json

ssl._create_default_https_context = ssl._create_unverified_context

def validate_directory(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run epigenetic clocks on simulation input files.")
    parser.add_argument("input_dir", type=validate_directory, help="Path to simulation files directory.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    input_dir = args.input_dir

    # You can now use `input_dir` in your existing logic, e.g.,
    # for file in os.listdir(input_dir): ...

def prepare_mData(beta_file):
    fH = open (beta_file,"r")
    lines = fH.readlines()
    mData = []
    tokens = lines[0].split("\t")
    header  = "ProbeID," + ",".join(tokens)

    for line in lines:
        tokens = line.split("\t")
        out = ",".join(tokens)
        mData.append(out)
    ofile = beta_file.replace(".txt",".prepared.txt")
    fOut = open (ofile,"w")
    for line in mData:
        fOut.write(line)
    print("DONE!!!")
    return

def read_file(prepared_file):
    data = pd.read_csv(prepared_file, sep=",")
    data.set_index("ID_REF",inplace=True)
    # create GeoData instance ...
    class geodata:
        def __init__(self, metadata=None,dnam=data, rna=None):
            self.dnam = data
            self.rna = rna
            self.metadata = metadata
    mdata = geodata(data)

    print("NOTICE: Done reading data from local file ...")
    return mdata
#
def run_clocks(mDNA):
    gallery = ModelGallery()
    #Note for warnings for missing data (default is imputation)...
    C = ['Horvathv1','Hannum','PhenoAge','Horvathv2','PEDBE','DunedinPACE']
    Horvathv1_results = gallery.get("Horvathv1").predict(mDNA)
    Hannum_results = gallery.get("Hannum").predict(mDNA)
    PhenoAge_results = gallery.get("PhenoAge").predict(mDNA)
    Horvathv2_results = gallery.get("Horvathv2").predict(mDNA)
    PEDBE_results = gallery.get("PEDBE").predict(mDNA)
    DunedinPACE_results = gallery.get("DunedinPACE").predict(mDNA)
    
    print("NOTICE: Done running the clocks ...")
    clocks_results = [Horvathv1_results,Hannum_results,PhenoAge_results,Horvathv2_results,PEDBE_results,DunedinPACE_results]
    combined_results = pd.concat(clocks_results, axis=1)
    combined_results.columns = C
    return combined_results
if __name__ == "__main__":
    #prepare DNAm files 
    path = sys.argv[1]#"path_to_files/"
    files = os.listdir(path)
    for file in files:
        if ".txt" in file:
            prepare_mData(path+file)
    # run clocks 
    files = os.listdir(path)
    for file in files:
        if ".prepared.txt" in file :
            print(file)
            mDNA = read_file(path+file)
            combined_results = run_clocks(mDNA)
            ofile = path+file.replace(".txt",".biolearn.csv")
            fOut = open (ofile,'w')
            fOut.write('Horvath,Hannum,Levine,skinHorvath,PedBE,DUNEDIN\n')# header that match Methyclock ...

            combined_results.to_csv(ofile,index=False)
            print("___________________________________________")
