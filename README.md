# Mine-N-Match (MNM): A NCBI Taxonomy and Sequence Data Mining Tool

This repository contains a set of Python functions designed to mine and process sequence data from the National Center for Biotechnology Information (NCBI) databases. It leverages the Biopython library for interacting with NCBI's Entrez API and other libraries such as `pygbif`, `pandas`, and `progressbar` for additional functionalities.

## Description

The core functionality of this package revolves around fetching sequence data (both nucleotide and protein) from NCBI based on user-defined queries and a list of species. It includes robust error handling and the ability to correct species names using the Global Biodiversity Information Facility (GBIF) backbone taxonomy, and can automatically record species synonyms.

**Key Features:**

*   **`ncbi_fetch_species()`:** Fetches species names from NCBI for a given higher-order taxonomic group and returns a dictionary containing taxonomic information.
*   **`ncbi_mine_seq_data()`:** Mines NCBI for sequence data based on a list of species and a query, saving results to CSV and FASTA files.

## Installation

1. **Clone the repository:**
   ```bash
    git clone https://github.com/VisualPhysiologyDB/mine_n_match

2. **Install dependencies:** [Make sure you are working in the repository directory from here-after]

   A. Create a Conda environment for Mine-N-Match (make sure you have [Conda](https://www.anaconda.com/) installed)
   ```bash
   conda create --name mnm_env python=3.11
   ```
   ### THEN
   ```bash
   conda activate mnm_env
   ```
   B. Use the 'requirements.txt' file to download base package dependencies for MNM
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1.  **Example:**
    Here is a minimal usage example for the `ncbi_mine_seq_data` function. This example uses the `species_list` variable, which is a list of species names, and `query` variable which is a search string for Entrez.

    ```python
    # Import functions for mining NCBI
    from mine_ncbi_functions import ncbi_fetch_species, ncbi_mine_seq_data 
    # Import json so we can load any existing 
    import json
    # Email for when we query NCBI
    email = "sethfrazer@ucsb.edu"  # Replace with your email
    ```

     ```python
      # Example usage:
      taxa = "Mammalia"
      rank = "class"
      limit = 500
      report_dir = 'taxonomy_data'
      out_file = "mammalia_taxonomy"
      
      species_data = ncbi_fetch_species(email, report_dir=report_dir, out=out_file, taxa=taxa, rank=rank, limit=limit, verbose=False)
     ```
  
     ```python
      term = f"(opsin[Title] OR rhodopsin[Title] OR OPN[Title] OR rh1[Title] OR rh2[Title] OR Rh1[Title] OR Rh2[Title]) NOT partial[Title] NOT voucher[All Fields] NOT kinase[All Fields] NOT kinase-like[All Fields] NOT similar[Title] NOT homolog[Title] NOT opsin-like[Title]"
  
      ncbi_query_df, query_report_dir = ncbi_mine_seq_data(email=email, job_label='ncbi_mammalia_opsins', out='ncbi_mammalia_opsins', species_list=species_list[140:143], taxa_dictionary=species_data)
     ```

## Notes

*   The functions include extensive error handling, especially for network-related issues when querying NCBI. If an error occurs, intermediate progress is often saved to files, allowing you to resume from where you left off.
*   The `ncbi_mine_seq_data()` function automatically creates a directory to store the results, labeled with the job name and current date/time.
*   Using an API key is highly recommended to increase the number of allowed requests per second and avoid potential issues with NCBI's rate limiting.
*   The `ncbi_fetch_species` function can take some time depending on the search parameters, especially for a large taxonomic group.

## Contributing

Contributions to this project are welcome. Please feel free to submit pull requests or open issues on the GitHub repository.

## License

All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies
