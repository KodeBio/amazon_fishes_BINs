# -*- coding: utf-8 -*-
# Code for processing BIN genetic data for biogeographic analysis

### ----------------------------
### 1. Import Required Packages
### ----------------------------
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
import pandas as pd
import numpy as np
import requests
import math
import re

print("Packages loaded successfully.")

### ----------------------------
### 2. Define Utility Functions
### ----------------------------

def which(condition):
    """Find indices"""
    return np.where(condition)[0]

def flattenlist(lista):
    """Flatten a list of lists into a single list."""
    return [item for sublist in lista for item in sublist] 

def gb_info(nucleotide_id, email):
    """
    Retrieve the country, publication title and voucher for a nucleotide sequence
    from GenBank.
    Args:
        nucleotide_id (str): GenBank accession ID.
        email (str): User's email for Entrez access.
    Returns:
        list: [country, paper title, voucher] or [None, None, None] if not found.
    """
    country_paper = [None, None, None]
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
    except Exception as e:
        print(f"Error fetching GenBank data for {nucleotide_id}: {e}")
        return country_paper

    genbank_text = handle.read()

    # Extract country information
    country_match = re.search(r'/country="([^"]+)"', genbank_text)
    if country_match:
        country = country_match.group(1).strip()
        country_paper[0] = re.sub(r'\s+', ' ', country)

    # Extract publication title
    title_match = re.search(r'TITLE\s+(.*\n(?:\s{12}.*\n)*)', genbank_text, re.MULTILINE)
    if title_match:
        country_paper[1] = re.sub(r'\s+', ' ', title_match.group(1).strip())

    # Extract specimen voucher information
    voucher_match = re.search(r'/specimen_voucher="([^"]+)"', genbank_text)
    if voucher_match:
        country_paper[2] = voucher_match.group(1).strip()
    
    return country_paper

print("Utility functions defined.")

### ----------------------------
### 3. Load and Process BIN Data
### ----------------------------

# Load Actinopterygii BIN data
file = "data/Actinopterygii_BINs.csv"
print(f"Loading BIN data from {file}...")
peces_bin = pd.read_csv(file, header=0, encoding='unicode_escape', low_memory=False)
print("BIN data loaded successfully.")

# Extract unique BIN URIs
BINs = list(set(peces_bin['bin_uri']))
print(f"Found {len(BINs)} unique BIN URIs.")

# Fill countries with BOLD information
print("Filling country information for BINs...")
for BIN in BINs:
    pos = which(peces_bin['bin_uri'] == BIN)
    df = peces_bin.iloc[pos]
    pais_in_df = df['country']

    if not any(pais_in_df.isna()):
        continue

    pos = pais_in_df.index.tolist()
    na_pos = which(pais_in_df.isna())
    indices = [pos[i] for i in na_pos]

    url = f"http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri={BIN}"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    paises = soup.select('.geoName')
    countries = [pais.text for pais in paises]

    new_country = np.array('~' + '/'.join(countries))
    new_country = np.tile(new_country, len(indices))

    peces_bin.loc[indices, 'country'] = new_country

print("Country information updated.")

# Replace blank spaces with NaN
peces_bin = peces_bin.applymap(lambda x: np.nan if isinstance(x, str) and x.strip() == "" else x)

# Save updated data
output_file = "data/Actinopterygii_BINs_country.csv"
peces_bin.to_csv(output_file, index=False, encoding='utf-8')
print(f"Updated BIN data saved to {output_file}.")

### ----------------------------
### 4. Extract Geographic Coordinates
### ----------------------------

print("Extracting geographic coordinates...")
peces_bin_py = peces_bin[~peces_bin['lat'].isna()]
paises = ["Bolivia", "Brazil", "Colombia", "Ecuador", "Peru", "Guyana", "Venezuela"]
sel_coun = peces_bin_py['country'].str.contains('|'.join(paises), case=False, na=False)
sel_coun = sel_coun | peces_bin_py['country'].isna()
peces_bin_py = peces_bin_py[sel_coun]

output_path = "data/amz_coords.txt"
peces_bin_py[['lat', 'lon']].to_csv(output_path, sep='\t', index=False)
print(f"Coordinates saved to {output_path}.")

### ----------------------------
### 5. Fetch GenBank Metadata
### ----------------------------

print("Fetching GenBank metadata...")
gb_code = list(set(peces_bin['genbank_accession'].dropna()))
papers_gb = pd.DataFrame({
    "gb_code": gb_code,
    "voucher": None,
    "gb_country": None,
    "paper": None
})

for i in range(papers_gb.shape[0]):
    code = papers_gb.loc[i, 'gb_code']
    
    info = gb_info(code, 'kiefer.bedoya@unmsm.edu.pe')
    
    papers_gb.loc[i, 'voucher'] = info[2]
    papers_gb.loc[i, 'gb_country'] = info[0]
    papers_gb.loc[i, 'paper'] = info[1]

    if (i + 1) % 10 == 0 or i == papers_gb.shape[0] - 1:
        print(f"Processed {i + 1}/{papers_gb.shape[0]} GenBank entries.")

papers_gb.sort_values(by='paper', inplace=True)
output_file = "data/genbank_data.csv"
papers_gb.to_csv(output_file, index=False, encoding='utf-8')
print(f"GenBank metadata saved to {output_file}.")

### ----------------------------
### 6. Fetch Habitat Information
### ----------------------------

print("Fetching habitat information...")
species = list(set(peces_bin['species_name'].dropna()))
species = [s.replace(" cf. ", " ") for s in species]
species = [' '.join(s.split()[:2]) for s in species]
species = list(set(species))

habitats = pd.DataFrame({
    "specie": species,
    "habitat": None
})

for i in range(habitats.shape[0]):
    url = f"https://www.fishbase.se/summary/{habitats.loc[i, 'specie'].replace(' ', '-')}.html"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    boxes = soup.find_all("div", class_="smallSpace")

    if not boxes:
        habitats.loc[i, 'habitat'] = "Unknown"
        continue

    mar_bool = re.search('marine', boxes[1].text, re.IGNORECASE)
    fresh_bool = re.search('freshwater', boxes[1].text, re.IGNORECASE)
    if mar_bool and fresh_bool:
        habitats.loc[i, 'habitat'] = "Marine and freshwater"
    elif fresh_bool:
        habitats.loc[i, 'habitat'] = "Freshwater"
    elif mar_bool:
        habitats.loc[i, 'habitat'] = "Marine"
    else:
        habitats.loc[i, 'habitat'] = "Unknown"

    if (i + 1) % 10 == 0 or i == habitats.shape[0] - 1:
        print(f"Processed {i + 1}/{habitats.shape[0]} habitat entries.")

output_file = "data/habitats.csv"
habitats.to_csv(output_file, index=False, encoding='utf-8')
print(f"Habitat information saved to {output_file}.")
