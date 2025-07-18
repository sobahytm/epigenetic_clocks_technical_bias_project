
import os,sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def validate_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path

def validate_directory(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path

def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze original and simulation output tables to compute delta age.")
    parser.add_argument("original_table", type=validate_file, help="CSV file containing the original predicted age table")
    parser.add_argument("simulation_tables_dir", type=validate_directory, help="Directory with simulated predicted age tables")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    original_table = args.original_table
    simulation_tables_dir = args.simulation_tables_dir



def read_table(file):
    df = pd.read_csv(file, sep=',', index_col='id')
    #df.set_index('id', inplace=True)
    return df
def analyze_datasets (original_df,sim_dfs):

  
    clocks_of_interest = ["Horvath", "Hannum", "Levine", "skinHorvath", "PedBE", "DUNEDIN"] 
    # (1) Calculate delta_age for each individual across 10 simulations
    delta_age_per_sim = []

    for sim_df in sim_dfs:
        delta_age = abs(sim_df[clocks_of_interest] - original_df[clocks_of_interest]) # using absolute delta age difference...
        delta_age_per_sim.append(delta_age)
    # Convert the list of delta_age DataFrames to a single DataFrame for averaging
    delta_age_avg_individuals = pd.concat(delta_age_per_sim).groupby(level=0).mean()

    delta_age_std_individuals = pd.concat(delta_age_per_sim).groupby(level=0).std()
    # (1.b) Generate a heatmap for the average delta age per individual
    #plt.figure(figsize=(10, 12))
    #sns.heatmap(delta_age_avg_individuals, cmap="coolwarm")
    #plt.title("American\nStress Plus Test\nAverage Delta Age \nper individual across simulations")
    #plt.show()

    # (2) Calculate the average delta age across individuals for each clock
    delta_age_avg_clocks = delta_age_avg_individuals.mean(axis=0) 
    delta_age_std_clocks = delta_age_std_individuals.mean(axis=0)
    # Convert it into a DataFrame for clarity
    delta_age_avg_clocks_df = pd.DataFrame(delta_age_avg_clocks, columns=['Average Delta Age'])
    delta_age_std_clocks_df = pd.DataFrame(delta_age_std_clocks, columns=['Average Delta Age'])
    # Display both tables
    print("Average Delta Age per Individual (Table 1):")
    print(delta_age_avg_individuals)

    print("\nAverage Delta Age per Epiclock (Table 2):")
    print(delta_age_avg_clocks_df)

    print("\nStd Delta Age per Epiclock (Table 3):")
    print(delta_age_std_clocks_df)

    delta_age_avg_individuals.to_csv("delta_age_avg_individuals.csv", index=True)
    delta_age_avg_clocks_df.to_csv("delta_age_avg_clocks.csv", index=True)
    delta_age_std_clocks_df.to_csv("delta_age_std_clocks.csv", index=True)
    return
#
# Adding a function to calculate percent deviation and integrating it into the analysis process

def calculate_percent_deviation(original_df, sim_dfs, clocks_of_interest):
    """
    Calculate percent deviation for each simulation compared to the original data.
    """
    percent_deviation_per_sim = []

    for sim_df in sim_dfs:
        percent_deviation = (
            (sim_df[clocks_of_interest] - original_df[clocks_of_interest]) / original_df[clocks_of_interest].replace(0, np.nan)) * 100
        percent_deviation_per_sim.append(percent_deviation)

    # Calculate the average percent deviation across simulations
    percent_deviation_avg_individuals = pd.concat(percent_deviation_per_sim).groupby(level=0).mean()

    # Calculate the average percent deviation for each clock
    percent_deviation_avg_clocks = percent_deviation_avg_individuals.mean(axis=0)

    # Convert the results into a DataFrame
    percent_deviation_avg_individuals.to_csv("percent_deviation_avg_individuals.csv", index=True)
    percent_deviation_avg_clocks_df = pd.DataFrame(percent_deviation_avg_clocks, columns=["Average Percent Deviation"])
    percent_deviation_avg_clocks_df.to_csv("percent_deviation_avg_clocks.csv", index=True)

    return percent_deviation_avg_individuals, percent_deviation_avg_clocks_df
if __name__ == "__main__":
    original_df = read_table(sys.argv[1])
    print('NOTICE: Done reading original file...')
    files = os.listdir(sys.argv[2])
    sim_dfs =[]
    for file in files:
        sim_df = read_table(sys.argv[2]+file)
        sim_dfs.append(sim_df)
        print('NOTICE: Done reading simulation file ...')

    analyze_datasets (original_df,sim_dfs)
    # Define clocks of interest (you can modify as needed)
    clocks_of_interest = ["Horvath", "Hannum", "Levine", "skinHorvath", "PedBE", "DUNEDIN"]
    
    # Call the new percent deviation function
    percent_deviation_avg_individuals, percent_deviation_avg_clocks_df = calculate_percent_deviation(
        original_df, sim_dfs, clocks_of_interest
    )

    # Print results for verification
    print("\nAverage Percent Deviation per Individual (Table 4):")
    print(percent_deviation_avg_individuals)

    print("\nAverage Percent Deviation per Clock (Table 5):")
    print(percent_deviation_avg_clocks_df)
    
