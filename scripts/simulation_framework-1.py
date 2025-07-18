
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
    parser = argparse.ArgumentParser(description="Run the simulation framework using mutation data and beta matrix.")
    parser.add_argument("mutation_file", type=validate_file, help="Mutation file: [population].common_mutations_in_CpG.with_zygosity.txt")
    parser.add_argument("intersected_data_dir", type=validate_directory, help="Directory containing intersected mutation data")
    parser.add_argument("beta_file", type=validate_file, help="DNA methylation beta matrix file (TSV format)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    # These will be passed to the function logic
    mutation_file = args.mutation_file
    intersected_data_dir = args.intersected_data_dir
    beta_file = args.beta_file


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
    cpgs =[]  # also create cpgs list for tracking ...
    files = os.listdir(path_intersected_data)
    for file in files:
        fH = open (path_intersected_data + file,'r')
        lines = fH.readlines()
        file_name = file.split('.')
        population = file_name[-4]
        if population in data:
            for line in lines:
                tokens = line.replace("\n","").split(",")
                CpG_mutations = tokens[-1].split(";")
                for CpG_mutation in CpG_mutations:
                    if CpG_mutation in data[population]:
                        data[population][CpG_mutation]["CpG"] = tokens[2]
                        data[population][CpG_mutation]["clock"].append(tokens[1])
                        cpgs.append(tokens[2])
    selected_per_cpg = {cpg: set() for cpg in cpgs}
    return data , selected_per_cpg
#
def read_DNAm_dataset(DNAm_dataset,data, zygosity_df,i,selected_per_cpg):
    # read data ...
    df = pd.read_csv(DNAm_dataset, index_col =0,sep = '\t')
    print('NOTICE: Done reading DNAm file...')
    # generate a list of the CpG overlap with a genetic variant and 
    # dictionary to map overlapping CpG to the genetic variant zygosity...
    # dictionary to map overlapping CpG to the Normalized_Probability...
    CpG_overlap_with_G = []
    CpG_to_Zygosity = {}
    CpG_to_Normalized_Probability = {}
    CpG_to_af  ={}
    for population in data:
        print(population)
        for variant in data[population]:
            CpG_overlap_with_G.append(data[population][variant]["CpG"])
            CpG_to_Zygosity[data[population][variant]["CpG"]] = data[population][variant]["zygosity"]
            CpG_to_af[data[population][variant]["CpG"]] = data[population][variant]["af"]
        #
        for variant in range(len(zygosity_df[population])):
            CpG_to_Normalized_Probability[data[population][zygosity_df[population][variant]]["CpG"]] = zygosity_df['Normalized_Probability'][variant]
    # create two lists ...
    cpgs = []
    n_probabilities = []
    items = CpG_to_Normalized_Probability.items()
    for item in items:
        cpgs.append(item[0]), n_probabilities.append(item[1])

    if len(list(set(CpG_overlap_with_G))) < len(CpG_overlap_with_G):
        print("WARNING: Please check the CpGs for duplicates!")
    CpG_overlap_with_G = list(set(CpG_overlap_with_G))
    track = {}
    for cpg in CpG_overlap_with_G:
        track[cpg] = []
    # create simulated dataset ...
    df_new = df.copy() # create a copy ...
    # Create a DataFrame to track which samples have been selected for each CpG
    tracking_df = pd.DataFrame(index=CpG_overlap_with_G, columns=df.columns[1:])  # Exclude 'ID_REF' column
    tracking_df[:] = False  # Initialize tracking DataFrame with False (no samples selected yet)
    
    # Loop over each CpG
    for cpg in CpG_overlap_with_G:
        available_samples = list(df.columns[1:])  # Get the sample column names
        # determine the selected sample size by AF ...
        sample_size =  {sample: set() for sample in df.columns}
        n = int(float(CpG_to_af[cpg]) * len(sample_size))
        if n < 5:#in case of low allele frequency ....
            n = 5 # minimum number of selected samples ...
        print('NOTICE: For {} {} {} will be selected'.format(cpg,n,float(CpG_to_af[cpg])))
        # Randomly select n samples from available samples for this CpG in this iteration
        selected_samples = np.random.choice(available_samples,n)
        # Mark the selected samples in the tracking DataFrame to avoid re-selection
        # check for duplication in the selected samples
        if len(selected_samples) != len(list(set(selected_samples))):
            selected_samples = list(set(selected_samples))
        tracking_df.loc[cpg, selected_samples] = True  
        # Remove selected samples from available samples for subsequent iterations
        available_samples = [sample for sample in available_samples if sample not in selected_samples]

        # Apply changes to the beta values based on zygosity
        for sample in selected_samples:
            # Check if CpG is present in the DataFrame
            if cpg in df_new.index:
                beta = df_new.at[cpg, sample]  # Get the current beta value
                track[cpg].append(sample)
                
                # Check zygosity and apply changes
                if beta >= 0:#0.3:  # Proceed only if beta value is greater than 0.3
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
    
    outfile = DNAm_dataset.replace('.txt','.simulated_iv5.{}.txt'.format(i))
    df_new.to_csv(outfile, sep='\t', index=True)
    # Verify the file creation
    print("NOTICE: DONE simulation {} ...".format(i))
    with open("track.simulation_i.{}.json".format(i), "w") as outfile: 
        json.dump(track, outfile)
    return

#

if __name__ == "__main__":
    data, zygosity_df = read_mutations_file(sys.argv[1])
    # add the clock information ... 
    updated_data, selected_per_cpg = process_intersected_data(sys.argv[2],data)
    print('NOTICE: Tracker is created ...')
    # generate simulated dataset ...
    for i in range (10): # set simulations iterations # 10 ...
        read_DNAm_dataset(sys.argv[3],data,zygosity_df, i,selected_per_cpg)
 

