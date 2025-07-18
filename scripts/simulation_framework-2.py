
import os, sys
import argparse
import pandas as pd
import numpy as np
import json

def validate_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path

def validate_directory(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the simulation framework (moderate/high intensity tests).")
    parser.add_argument("mutation_file", type=validate_file, help="Mutation file: [population].common_mutations_in_CpG.with_zygosity.txt")
    parser.add_argument("intersected_data_dir", type=validate_directory, help="Directory containing intersected mutation data")
    parser.add_argument("beta_file", type=validate_file, help="DNA methylation beta matrix file (TSV format)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    # These variables can now be used in the original function logic
    mutation_file = args.mutation_file
    intersected_data_dir = args.intersected_data_dir
    beta_file = args.beta_file

    # You may now call the original functions using these variables
    # Example:
    # mutations_data = read_mutations_file(mutation_file)

def read_mutations_file(mutation_file):
    file_name = mutation_file.split("/")
    file_name2 = file_name[-1].split(".")
    population = file_name2[0]
    fH = open (mutation_file,"r")
    lines = fH.readlines()
    del lines[0]
    data = {}
    data[population] = {}
    for line in lines:
        tokens = line.replace("\n","").split("\t")
        data[population][tokens[0]] = {"af":0,"zygosity":"na","CpG":"na","clock":[]}
        data[population][tokens[0]]["af"] = float(tokens[1])
        data[population][tokens[0]]["zygosity"] = tokens[2]
    # create adjusted probabilites ...
    zygosity_df = pd.read_csv(mutation_file,sep = '\t')

    # Define zygosity weight mapping
    zygosity_weight = {
        'het': 1,    # Heterozygous variants
        'homo': 2    # Homozygous variants (if present)
    }

    # Apply zygosity weights to the allele frequencies
    zygosity_df['Adjusted_AF'] = zygosity_df.apply(lambda row: row['AF'] * zygosity_weight.get(row['zyogosity_probability'], 1), axis=1)

    # Normalize the adjusted allele frequencies
    total_adjusted_AF = zygosity_df['Adjusted_AF'].sum()
    zygosity_df['Normalized_Probability'] = zygosity_df['Adjusted_AF'] / total_adjusted_AF



    return data, zygosity_df
#
def process_intersected_data(path_intersected_data, data):
    files = os.listdir(path_intersected_data)
    for file in files:
        fH = open (path_intersected_data + file,'r')
        lines = fH.readlines()
        file_name = file.split('.')
        population = file_name[-2]# use file_name[-4], if common SNPs 
        if population in data:
            for line in lines:
                tokens = line.replace("\n","").split(",")
                CpG_mutations = tokens[-1].split(";")
                for CpG_mutation in CpG_mutations:
                    if CpG_mutation in data[population]:
                        data[population][CpG_mutation]["CpG"] = tokens[2].replace('*','')
                        data[population][CpG_mutation]["clock"].append(tokens[1])
    return data
#
def read_DNAm_dataset(DNAm_dataset,data, zygosity_df,i,selected_per_sample):
    # read data ...
    df = pd.read_csv(DNAm_dataset, index_col =0,sep = '\t')
    print('NOTICE: Done reading DNAm file...')
    # generate a list of the CpG overlap with a genetic variant and 
    # dictionary to map overlapping CpG to the genetic variant zygosity...
    # dictionary to map overlapping CpG to the Normalized_Probability...
    CpG_overlap_with_G = []
    CpG_to_Zygosity = {}
    CpG_to_Normalized_Probability = {}
    for population in data:
        print(population)
        for variant in data[population]:
            CpG_overlap_with_G.append(data[population][variant]["CpG"])
            CpG_to_Zygosity[data[population][variant]["CpG"]] = data[population][variant]["zygosity"]
        #
        for variant in range(len(zygosity_df[population])):
            CpG_to_Normalized_Probability[data[population][zygosity_df[population][variant]]["CpG"]] = zygosity_df['Normalized_Probability'][variant]
    # create  lists ...
    n_probabilities = []
    items = CpG_to_Normalized_Probability.items()
    for item in items:
        CpG_overlap_with_G.append(item[0]), n_probabilities.append(item[1])

    if len(list(set(CpG_overlap_with_G))) < len(CpG_overlap_with_G):
        print("WARNING: Please check the CpGs for duplicates!")
        CpG_overlap_with_G = list(set(CpG_overlap_with_G))
    # create simulated dataset ...
    # create tracking dictionary ...
    track = {}
    for cpg in CpG_overlap_with_G:
        track[cpg] = []
    df_new = df.copy() # create a copy ...
    not_found = []#Oct27 ...
    # Iterate through each sample 
    for sample in df_new.columns:
        # Get the CpGs that haven't been selected yet for this sample
        available_CpGs = [cpg for cpg in CpG_overlap_with_G if cpg not in selected_per_sample[sample]]
       # set samples selection size ...
        n = 0.07 * 656 # Adjust per DNAm dataset (use mean of AFs * sample size)
        if len(available_CpGs) >= n:
            n = int(n)
            selected_CpGs = np.random.choice(available_CpGs, n, replace=False)
        else:
            selected_CpGs = available_CpGs
        # Update the set of previously selected CpGs for the current sample
        selected_per_sample[sample].update(selected_CpGs)

        # Apply changes to the beta values based on zygosity
        for cpg in selected_CpGs:
            # Check if CpG is present in the DataFrame
            if cpg in df_new.index:# if cpg not here it is not in this array ...
                beta = df_new.at[cpg, sample]  # Get the current beta value
                track[cpg].append(sample)
                
                # Check zygosity and apply changes
                if beta >= 0:  # Proceed only if beta value is greater than 0.3
                    if CpG_to_Zygosity.get(cpg) == "het":
                        # Heterozygous case: reduce by a percentage between 0% and 50%
                        percentage_to_change = np.random.uniform(0.0, 0.5)
                    else:  # Homozygous case: reduce by a percentage between 60% and 100%
                        percentage_to_change = np.random.uniform(0.6, 1.0)
                    
                    change = beta * percentage_to_change
                    new_beta = beta - change
                    
                    # If new beta is negative, print a warning and exit
                    if new_beta < 0:
                        print(f"WARNING: Beta value is negative for {cpg} in sample {sample}: {new_beta}")
                        exit(0)
                    
                    # Update the beta value in the DataFrame
                    df_new.at[cpg, sample] = new_beta
            else:# Oct27
                missing_cpgs_from_450k = ['cg06094762','cg08724636','cg10959651','cg11620135','cg14361627','cg17238334','cg18769120','cg20674577','cg21944491','cg22029879','cg22512531','cg23091758','cg26311454','cg26665419']
                if cpg in missing_cpgs_from_450k:
                    not_found.append(cpg)
    not_found = list(set(not_found))# Oct27
    print(not_found)# Oct27
    outfile = DNAm_dataset.replace('.txt','.simulated_ii.{}.txt'.format(i))
    df_new.to_csv(outfile, sep='\t', index=True)
    # Verify the file creation
    print("NOTICE: DONE simulation {} ...".format(i))
    with open("track.simulation.{}.json".format(i), "w") as outfile: 
        json.dump(track, outfile)
    return
#
def develope_tracker(DNAm_file):
    df = pd.read_csv(DNAm_file, index_col =0,sep = '\t')
    cpgs = list(df.index)  # List of CpGs, assuming CpGs are the row indices of the DataFrame
    # Initialize a dictionary to track previously selected CpGs for each sample
    selected_per_sample = {sample: set() for sample in df.columns}
    print('NOTICE: Created tracker ...')
    return selected_per_sample
#

if __name__ == "__main__":
    data, zygosity_df = read_mutations_file(sys.argv[1]);print()
    # add the clock information ... 
    updated_data = process_intersected_data(sys.argv[2],data)
    # create tracker  ...
    selected_per_sample = develope_tracker(sys.argv[3])
    # generate simulated dataset ...
    for i in range (3): # 10 simulations ...
        read_DNAm_dataset(sys.argv[3],data,zygosity_df, i,selected_per_sample)

