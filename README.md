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
│   ├── enrichment_analysis_clocks.py 
|   |── enrichment_analysis_populations.py
|   |── clocks_weights_based_analysis.py
│   ├── simulation_framework-1.py
│   ├── simulation_framework-2.py
│   ├── run_clocks.py
│   └── post-simulation_analysis.py
├── data/
│   ├── summary_CpG_all_snps.txt
|   |── gnomAD_overlap_table.snps.csv
|   |── clock_coefficient.txt
|   |── intersected_data/
|   |── common_snps/
|   |── all_snps/
├── results/             # Optional (example outputs)
├── LICENSE
├── .gitignore
└── README.md
```

---

## 🚀 Usage

Each script is independently executable from the command line using Python 3.8+.


### 1. Enrichment Analysis: Clocks
```bash
python scripts/Enrichment_analysis_clocks.py \
    path/to/summary_CpG_all_snps.txt
```

### 2. Enrichment Analysis: Populations
```bash
python scripts/enrichment_analysis_populations.py \
    path/to/gnomAD_overlap_table.snps.csv
```

### 3. CpG Weight-Based Enrichment
```bash
python scripts/Clocks_weights_based_analysis.py \
    path/to/clock_coefficient.txt \
    path/to/intersected_data/
```

### 4. Simulation Test I
```bash
python scripts/simulation_framework-1.py \
    path/to/common_snps/[population].with_zygosity.txt \
    path/to/intersected_data/ \
    path/to/beta_matrix.txt
```

### 5. Simulation Test II
```bash
python scripts/simulation_framework-2.py \
    path/to/all_snps/[population].with_zygosity.txt \
    path/to/intersected_data/ \
    path/to/beta_matrix.txt
```

### 6. Run Epigenetic Clock Models
```bash
python scripts/run_clocks.py \
    path/to/DNA_methylation_directory/
```

### 7. Post-Simulation Analysis
```bash
python scripts/post-simulation_analysis.py \
    path/to/original_clock_output.csv \
    path/to/simulated_outputs_directory/
```

---


## ⚖️ License

This project will be released under the [MIT License](LICENSE).

---

## 📚 Citation

The associated manuscript is currently in preparation and will be added upon publication.
