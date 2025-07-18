
"""
Epigenetic Clock Enrichment Analysis
------------------------------------
This script computes the enrichment of genetic mutations in CpG sites used by epigenetic clocks across populations.
It includes:
- Fisher's exact test per clock vs. others per population
- FDR-adjusted p-values (Benjamini-Hochberg)
- Relative Risk (RR) calculation
- Heatmap generation for both metrics
- CSV export of p-values

"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def load_data(filepath):
    df = pd.read_csv(filepath)
    exclude = ['Monika', 'Carola', 'Maria']
    df = df[~df['clock'].isin(exclude)].reset_index(drop=True)
    return df

def compute_fisher_with_fdr(df, populations):
    population_columns = {
        pop: (f"{pop}_site_with_mutations", f"{pop}_site_without_mutations") for pop in populations
    }
    clock_names = df['clock'].tolist()
    raw_p = []
    fisher_results = {pop: [] for pop in populations}

    for pop in populations:
        mut_col, nomut_col = population_columns[pop]
        for i in range(len(clock_names)):
            this_mut = df.loc[i, mut_col]
            this_nomut = df.loc[i, nomut_col]
            others_mut = df[mut_col].sum() - this_mut
            others_nomut = df[nomut_col].sum() - this_nomut
            table = np.array([[this_mut, this_nomut], [others_mut, others_nomut]])
            _, p = fisher_exact(table, alternative='two-sided')
            raw_p.append(p)
            fisher_results[pop].append(p)

    # Adjust p-values
    adj_p = multipletests(raw_p, method='fdr_bh')[1]

    # Build result table
    table = []
    i = 0
    for pop in populations:
        for clock in clock_names:
            table.append([pop, clock, raw_p[i], adj_p[i]])
            i += 1
    pval_df = pd.DataFrame(table, columns=["Population", "Clock", "Raw_p", "Adjusted_p"])

    # Matrix for heatmap
    matrix = []
    index = 0
    for pop in populations:
        row = []
        for _ in clock_names:
            val = adj_p[index]
            row.append(-np.log10(val) if val > 0 else 0)
            index += 1
        matrix.append(row)
    matrix_df = pd.DataFrame(matrix, index=populations, columns=clock_names)
    matrix_df.loc['average'] = matrix_df.mean(axis=0)

    return pval_df, matrix_df

def compute_relative_risk(df, populations):
    rr_matrix = []
    clock_names = df['clock'].tolist()
    population_columns = {
        pop: (f"{pop}_site_with_mutations", f"{pop}_site_without_mutations") for pop in populations
    }

    for pop in populations:
        mut_col, nomut_col = population_columns[pop]
        row = []
        for i in range(len(clock_names)):
            a = df.loc[i, mut_col]
            b = df.loc[i, nomut_col]
            a_total = a + b

            others_mut = df[mut_col].sum() - a
            others_total = others_mut + df[nomut_col].sum() - b

            clock_risk = a / a_total if a_total > 0 else 0
            others_risk = others_mut / others_total if others_total > 0 else 0
            rr = clock_risk / others_risk if others_risk > 0 else np.nan
            row.append(rr)
        rr_matrix.append(row)
    rr_df = pd.DataFrame(rr_matrix, index=populations, columns=clock_names)
    rr_df.loc['average'] = rr_df.mean(axis=0)
    return rr_df

def plot_heatmap(dataframe, title, center_val, outname):
    plt.figure(figsize=(12, 6))
    sns.heatmap(
        dataframe.T,
        cmap='coolwarm',
        linewidths=0.8,
        annot=True,
        fmt=".2f",
        center=center_val,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"orientation": "horizontal", "shrink": 0.8}
    )
    plt.xlabel("Population")
    plt.ylabel("Epigenetic Clock")
    plt.xticks(rotation=45, ha='right')
    plt.title(title)
    plt.tight_layout()
    #plt.savefig(outname, dpi=300)
    #plt.close()
    plt.show()
def main():
    parser = argparse.ArgumentParser(description="Epigenetic Clocks Mutation Enrichment Analysis")
    parser.add_argument("-i", "--input", required=True, help="Input summary file")
    parser.add_argument("-o", "--outdir", required=False, default="results", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_data(args.input)
    populations = ['afr', 'sas', 'amr', 'eas', 'fin', 'nfe', 'asj', 'ami', 'mid']

    # Fisher test + FDR
    pval_df, log_fdr_df = compute_fisher_with_fdr(df, populations)
    pval_df.to_csv(os.path.join(args.outdir, "fisher_fdr_pvalues.csv"), index=False)
    plot_heatmap(log_fdr_df, "-log10(FDR-adjusted p-values)", center_val=0, outname=os.path.join(args.outdir, "fdr_heatmap.png"))

    # Relative risk
    rr_df = compute_relative_risk(df, populations)
    rr_df.to_csv(os.path.join(args.outdir, "relative_risk_scores.csv"))
    plot_heatmap(rr_df, "Relative Risk", center_val=1, outname=os.path.join(args.outdir, "rr_heatmap.png"))

if __name__ == "__main__":
    main()
