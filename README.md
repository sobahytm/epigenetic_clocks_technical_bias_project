# Investigating the Influence of SNP-CpG Overlaps on Epigenetic Clock Performance in Diverse Populations

Epigenetic clocks, built from DNA methylation patterns, are among the most accurate biomarkers of aging. These models are typically trained using array-based platforms that measure methylation at hundreds of thousands of CpG sites across the genome. However, underlying genetic variation, particularly single-nucleotide polymorphisms (SNPs) overlapping CpG sites, can introduce technical noise into methylation measurements. Because allele frequencies vary across populations, such SNP-CpG overlaps may differentially affect model predictions.

To address this, we systematically evaluated the impact of SNP-CpG site overlaps on epigenetic clock performance in diverse populations using three complementary approaches:
- **Enrichment analysis**
- **CpG site weight evaluation**
- **Simulation-based analysis**

---

## 📁 Repository Structure

```
epigenetic-clocks-simulation/
├── scripts/
│   ├── Clocks_weights_based_analysis_polished.py
│   ├── Simulation_framework_polished.py
│   ├── Simulation_framework-2_polished.py
│   ├── Run_clocks_polished.py
│   └── Post-simulation_analysis_polished.py
├── data/                # Optional (sample input files or instructions)
├── results/             # Optional (example outputs)
├── docs/                # Optional (extended documentation or figures)
├── LICENSE
├── .gitignore
└── README.md
```

---

## 🚀 Usage

Each script is independently executable from the command line using Python 3.8+.

### 1. Enrichment Analysis: Populations
```bash
python scripts/Enrichment_analysis.Populations_polished.py \
    path/to/variants_overlap_table.snps.csv
```

### 2. Enrichment Analysis: Clocks
```bash
python scripts/Enrichment_Clocks_Enhanced.py \
    path/to/summary_CpG_all_mutations.txt
```


### 3. CpG Weight-Based Enrichment
```bash
python scripts/Clocks_weights_based_analysis_polished.py \
    path/to/clocks_coefficient.txt \
    path/to/mutation_directory/
```

### 4. Simulation Test I
```bash
python scripts/Simulation_framework_polished.py \
    path/to/population.common_mutations_in_CpG.with_zygosity.txt \
    path/to/intersected_data_directory/ \
    path/to/beta_matrix.tsv
```

### 5. Simulation Test II
```bash
python scripts/Simulation_framework-2_polished.py \
    path/to/population.common_mutations_in_CpG.with_zygosity.txt \
    path/to/intersected_data_directory/ \
    path/to/beta_matrix.tsv
```

### 6. Run Epigenetic Clock Models
```bash
python scripts/Run_clocks_polished.py \
    path/to/simulation_results_directory/
```

### 7. Post-Simulation Analysis
```bash
python scripts/Post-simulation_analysis_polished.py \
    path/to/original_clock_output.csv \
    path/to/simulated_outputs_directory/
```

---



## ⚖️ License

This project will be released under the [MIT License](LICENSE).

---

## 📚 Citation

The associated manuscript is currently in preparation and will be added upon publication.
