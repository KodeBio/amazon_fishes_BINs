# **Amazon_Fishes_BINs**

**Repository for Manuscript:**
*Advances in Biogeographic Analysis of Genetic Data Using Database-Driven Metrics: A Case Study on the Amazon Basin*

This repository contains the data, scripts, and procedures necessary to reproduce the analyses presented in the manuscript. The workflow integrates geographic and genetic data using database-driven methods to advance biogeographic research within the Amazon Basin.

---

## **Workflow Overview**
The following steps describe the complete workflow, from data acquisition to final analysis:

### **1. Geographic Data Acquisition**
Download the **HydroBasin** geographic dataset from HydroSHEDS:
- **Website**: [HydroSHEDS - HydroBasins](https://www.hydrosheds.org/products/hydrobasins)
- **Data Required**: Standard datatype for South America.
- **Output**: `hybas_sa_lev04_v1c.shp`

This shapefile contains the hydrological basins of South America, essential for spatial analysis of genetic data.

---

### **2. Genetic Data Acquisition**
Download genetic data from the BIN database:
- **Website**: [BIN Database](https://is.gd/7Ea37q)
- **Target Taxonomic Group**: *Actinopterygii* (ray-finned fishes), worldwide.
- Convert the downloaded data into a CSV file format and save it as:
  - `Actinopterygii_BINs.csv`

This dataset provides the genetic BIN records for subsequent analysis.

---

### **3. Extracting Genetic and Geographic Information**
Run the script `Getting_information.py` using `Actinopterygii_BINs.csv` as input. This script processes the genetic dataset to generate the following outputs:
- **`Actinopterygii_BINs_country.csv`**: Contains BIN records annotated with country information directly obtained from the webpage of BOLD.
- **`amz_coords.txt`**: Contains the geographic coordinates of BIN sequences.
- **`genbank_data.csv`**: Includes metadata from GenBank for sequences in the BIN database.
- **`habitats.csv`**: Lists habitat types for each BIN sequence.

---

### **4. Assigning Basin Information**
Run the script `Location_by_coords.py`:
- **Input Files**:
  - `amz_coords.txt` (generated in Step 3)
  - `hybas_sa_lev04_v1c.shp` (downloaded in Step 1)
- **Output File**:
  - `amz_basin.txt`: Provides basin-level information for each geographic coordinate.

---

### **5. Filtering and Data Preparation**
Run the script `Filtering.R` using all files created in previous steps:
- **Input Files**:
  - `Actinopterygii_BINs_country.csv`
  - `amz_basin.txt`
  - `genbank_data.csv`
  - `habitats.csv`
- **Output Files**:
  - `amz_BOLD_bins.csv`: Filtered BIN records for Amazon Basin.
  - `amz_BOLD_seqs.fasta`: FASTA file containing sequences for further alignment.

---

### **6. Sequence Alignment and Manual Curation**
Align the sequences in `amz_BOLD_seqs.fasta`:
1. Use **BioEdit** or another tool to manually inspect and remove poorly aligned sequences.
2. Save the curated alignment as:
   - `amz_BOLD_alig.fasta`

---

### **7. Biogeographic Analysis**
Run the script `Analysis.R`:
- **Input Files**:
  - `amz_BOLD_bins.csv`
  - `amz_BOLD_alig.fasta`
- **Outputs**:
  - Analytical results, including visualizations and metrics.

The script uses various analytical functions stored in `utils.R` to perform the core biogeographic and genetic analyses.

---
