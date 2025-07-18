
import sys
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt

def compute_relative_risk(a, b, c, d):
    if b == 0 or d == 0:
        return np.nan
    pop_risk = a / b
    if c == 0:
        other_risk = (c + 0.5) / d  # continuity correction
    else:
        other_risk = c / d
    return pop_risk / other_risk if other_risk > 0 else np.nan

def compute_fisher_pvalue(a, b, c, d):
    contingency_table = [[a, b - a], [c, d - c]]
    _, p_value = fisher_exact(contingency_table, alternative='two-sided')
    return p_value

def analyze_enrichment(input_file, output_csv, heatmap_path=None):
    df = pd.read_csv(input_file)

    rr_list = []
    pval_list = []

    for _, row in df.iterrows():
        a = row["a (SNPs in pop overlap CpGs)"]
        b = row["b (All SNPs in pop)"]
        c = row["c (Other pop SNPs overlap CpGs)"]
        d = row["d (All SNPs in other pops)"]

        rr = compute_relative_risk(a, b, c, d)
        pval = compute_fisher_pvalue(a, b, c, d)

        rr_list.append(rr)
        pval_list.append(pval)

    df["RR"] = rr_list
    df["p_value"] = pval_list

    # Compute FDR
    df["FDR"] = multipletests(df["p_value"], method='fdr_bh')[1]

    # Save result table
    #df.to_csv(output_csv, index=False)

    # Optional heatmap
    if heatmap_path:
        heatmap_data = df.pivot(index="Population", columns="Clock", values="RR")
        heatmap_data["Average"] = heatmap_data.mean(axis=1)
        plt.figure(figsize=(12, 6))
        sns.heatmap(
            heatmap_data,
            cmap="coolwarm",
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            center=1,
            cbar_kws={"orientation": "horizontal", "shrink": 0.8}
        )
        plt.title("Relative Risk (RR) Heatmap with Population Averages")
        plt.xlabel("Clock")
        plt.ylabel("Population")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        #plt.savefig(heatmap_path)
        #plt.close()
        plt.show()
if __name__ == "__main__":
    analyze_enrichment(
        input_file=sys.argv[1],#"variants_overlap_table.snps.csv",
        output_csv="enrichment_results_with_rr_fdr.csv",
        heatmap_path="rr_heatmap_corrected.png"
    )
