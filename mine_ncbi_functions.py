import os
import re
import datetime
import time
import json
import copy
import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
import progressbar

api_key = "1efb120056e1cea873ba8d85d6692abd5d09"

# create a function that return the collection of all NCBI fetch result
def ncbi_fetch(email, term, ncbi_db="nuccore", rettype="gb", format="genbank"):
    
  """Fetches sequences from NCBI's databases using Entrez.

    This function queries an NCBI database (default: "nuccore") using a provided search term and 
    retrieves all corresponding sequence records. It uses the provided email for identification 
    and an API key for increased request limits.

    Args:
        email (str): Your email address. Required by NCBI for Entrez queries.
        term (str): The search term to query the NCBI database.
        ncbi_db (str, optional): The NCBI database to search. Defaults to "nuccore".
        rettype (str, optional): The retrieval type for efetch. Defaults to "gb".
        format (str, optional): The format of the sequence records to retrieve. Defaults to "genbank".

    Returns:
        list: A list of Biopython SeqRecord objects representing the fetched sequence records.

    Notes:
        - An API key should be set using `api_key = 'your_api_key'` in the main code before calling this function for increased requests per second.
        - The function automatically pauses for 0.25 seconds between fetches to avoid overloading 
          the NCBI servers, adhering to the rate limit even with an API key.
        - The function returns all results, even if the number of results exceeds the default `retmax` limit of Entrez.esearch.
  """
  Entrez.email = email    # Always tell NCBI who you are
  handle = Entrez.esearch(db=ncbi_db,
                        term=term, 
                        api_key = api_key) # using api key allows 10 fetch per second
  record = Entrez.read(handle)
  #print(record)
  full_res = int(record["Count"])
  
  handle_full = Entrez.esearch(db=ncbi_db,
                        term=term, 
                        api_key = api_key,
                        retmax = full_res + 1)
  record_full = Entrez.read(handle_full)
  #print(record_full)
  q_list = record_full["IdList"]

# create a list for all the entries fetched from NCBI
  SeqRecords =[]
  for x in q_list:
    # fetch result from previous search
    fet = Entrez.efetch(db=ncbi_db, id=x, rettype=rettype)
    seq = SeqIO.read(fet, format)
    #print(seq)
    fet.close()
    
    SeqRecords.append(seq)

    time.sleep(0.25)  # Rate limiting
  
  return SeqRecords

from pygbif import species
def correct_species_name(species_name):
    
    """Attempts to find a corrected species name and its associated higher taxa using the GBIF backbone taxonomy.

    This function uses the `pygbif` library to query the Global Biodiversity Information Facility (GBIF) 
    backbone taxonomy. It attempts to find a match for the input `species_name` and returns the 
    corrected species name (if different) along with a higher taxonomic group (family, order, or genus).
    The higher taxonomic group can be used to find related taxa from NCBI if a species name is not recognized.

    Args:
        species_name (str): The species name to check.

    Returns:
        tuple: A tuple containing two strings:
            - The corrected species name (if found), or "No match found" otherwise.
            - The associated higher taxa (family, order, or genus) of the corrected species name (if found), or "No match found" otherwise.

    Notes:
        - The function retries the GBIF query up to 10 times if it fails initially.
        - The function prioritizes returning the family, then the order, then the genus if a species is found.
        - If no match is found in the GBIF backbone, or if there's an error during the query, 
          it returns ("No match found", "No match found").
    """
    queried = False
    try:
        for x in range(10):
            if queried == False:
                result = species.name_backbone(name=species_name, rank='species', verbose=True)
                queried = True
    except:
        pass
    
    try:
        if 'species' in result and 'family' in result:
            return result['species'], result['family']
        elif ('species' in result) and ('order' in result):
            return result['species'], result['order']
        elif ('species' in result):
            return result['species'], result['genus']
        else:
            return "No match found", "No match found"
    except:
        return "No match found", "No match found"



def get_species_taxonomy(species_name, email, record_alt=False, use_higher_taxa=False, higher_taxa=None):
    
    """Fetches taxonomic information, including synonyms, for a given species or higher taxon from NCBI Taxonomy.

    This function retrieves taxonomic details from the NCBI Taxonomy database using the Biopython Entrez module.
    It can retrieve synonyms for a species or information about a higher taxon if `use_higher_taxa` is set to True.
    It also attempts to determine the Phylum, Class, and Order of the species or taxon.

    Args:
        species_name (str): The scientific name of the species to look up.
        email (str): Your email address, required by NCBI Entrez.
        record_alt (bool, optional): If True, includes the input `species_name` in the synonyms list. Defaults to False.
        use_higher_taxa (bool, optional): If True, searches for a higher taxon instead of a species. Defaults to False.
        higher_taxa (str, optional): The name of the higher taxon to search for when `use_higher_taxa` is True.

    Returns:
        dict: A dictionary containing the taxonomic information, including:
            - Synonyms (list): A list of synonymous names for the species.
            - Phylum (str): The phylum of the species/taxon.
            - Class (str): The class of the species/taxon.
            - Order (str): The order of the species/taxon.
            - If a rank is not found, it will be assigned "Unknown".
    """

    Entrez.email = email    # Always tell NCBI who you are
    queried = False
    synonyms = []

    for x in range(10):
        if queried == False:
            try:
                if use_higher_taxa == False:
                    # Search for the species
                    handle = Entrez.esearch(db="taxonomy", term=species_name, api_key = api_key)
                else:
                    handle = Entrez.esearch(db="taxonomy", term=higher_taxa, api_key = api_key)

                record = Entrez.read(handle)
                taxonomy = {}
            
                tax_id = record["IdList"][0]
                # Fetch the taxonomy record
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                record = Entrez.read(handle)[0]
                #if the code makes it this far, then there was a successful query
                queried = True
                #print(record)
                # Extract synonyms
                if ("OtherNames" in record) and (use_higher_taxa == False):
                    for name in record["OtherNames"]["Synonym"]:
                        if name not in synonyms:
                            if '(' in name:
                                if ',' in name:
                                    rename = name.split('(')[0].strip()
                                    if (rename not in synonyms) and (rename != species_name): 
                                        synonyms.append(rename)
                                else:
                                    rename = name.split('(')[1].strip()
                                    rename = rename.split(')')[0].strip() + ' ' + rename.split(')')[1].strip()
                                    if (rename not in synonyms) and (rename != species_name):
                                        synonyms.append(rename)
                                    rename2 = name.split('(')[0].strip() + rename.split(' ')[1]
                                    if (rename2 not in synonyms) and (rename2 != species_name):
                                        synonyms.append(rename2)
                            else:
                                if ',' in name:
                                    rename = name.split(',')[0]
                                    rename = rename.split(' ')[0] + ' ' + rename.split(' ')[1]
                                    if (rename not in synonyms) and (rename != species_name):
                                        synonyms.append(rename)     
                                else:
                                    synonyms.append(name)
                        else:
                            pass
    
                for lineage in record["LineageEx"]:
                    rank = lineage["Rank"]
                    if rank == "phylum": 
                        taxonomy["Phylum"] = lineage["ScientificName"]
                    elif rank == "class":
                        taxonomy["Class"] = lineage["ScientificName"]
                    elif rank == "Order":
                        taxonomy["Order"] = lineage["ScientificName"]
            except:
                pass
            
    if record_alt == True:
        synonyms.append(species_name)
        
    # Asign species synonyms  
    taxonomy['Synonyms'] = synonyms
    
    # Ensure all desired ranks are present
    for rank in ["Phylum", "Order", "Class"]:
        if rank not in taxonomy:
            taxonomy[rank] = "Unknown"
            
    handle.close()
    
    return taxonomy

def get_sp_taxon_dict(species_list, email, taxon_file, sp_taxon_dict={}):
    
    """Builds a dictionary of taxonomic information for a list of species.

    This function takes a list of species names and constructs a dictionary where keys are 
    species names and values are dictionaries containing taxonomic information 
    (e.g., Phylum, Subphylum, Class, Order, Synonyms) for each species. It uses the 
    `get_species_taxonomy` function to retrieve the information from NCBI, and can also
    leverage the `correct_species_name` function to handle potential misspellings or outdated names.

    Args:
        species_list (list): A list of species names (strings).
        email (str): Your email address, required for NCBI Entrez queries.
        taxon_file (str): The file path to save the resulting taxonomy dictionary (JSON format). This is used in error handling to save progress.
        sp_taxon_dict (dict, optional): An existing taxonomy dictionary to update. Defaults to an empty dictionary.

    Returns:
        dict: A dictionary where keys are species names and values are dictionaries 
              containing taxonomic information for each species.

    Raises:
        Exception: If there is an issue querying NCBI for taxonomic information. In this case the
                   current progress is saved to the taxon_file, and the user is prompted to restart.
                   This is implemented to avoid issues with querying NCBI's servers.

    Notes:
        - If a species is already present in `sp_taxon_dict`, its information is not fetched again,
        unless the existing entry is missing crucial taxonomic levels, in which case it will try to look
        up the species again using the `correct_species_name` function if possible.
        - The function attempts to correct species names using the `correct_species_name` function if the initial
          taxonomy lookup fails or returns incomplete information (missing Phylum, Order, and Class).
        - The `correct_species_name` function tries to find a corrected name or higher taxonomic
          group from the Global Biodiversity Information Facility (GBIF).
        - Intermediate results are saved to `taxon_file` in case of errors during the process. This allows
          you to resume the process from where it left off if an error occurs, rather than starting from scratch.
    """
    all_keys = sp_taxon_dict.keys()

    for species in species_list:
        #if (species in all_keys) and (sp_taxon_dict[species]['Phylum']!="Unknown"):
        if (species in all_keys):
            pass
        else: 
            try:
                taxonomy = get_species_taxonomy(species, email)
                #print(f"{species} synonyms: {synonyms}")
                sp_taxon_dict[species] = taxonomy
                if (len(sp_taxon_dict[species]['Synonyms']) == 0) and (sp_taxon_dict[species]['Phylum']=="Unknown") and (sp_taxon_dict[species]['Order']=="Unknown") and (sp_taxon_dict[species]['Class']=="Unknown"):
                    try:
                        corrected_name, higher_taxa_name = correct_species_name(species)
                        if (corrected_name != "No match found") and (corrected_name != species):
                            taxonomy = get_species_taxonomy(corrected_name, email, record_alt=True)
                            sp_taxon_dict[species] = taxonomy

                            if (sp_taxon_dict[species]['Phylum']=="Unknown") and (sp_taxon_dict[species]['Order']=="Unknown") and (sp_taxon_dict[species]['Class']=="Unknown"):
                                try:
                                    taxonomy = get_species_taxonomy(corrected_name, email, record_alt=True, use_higher_taxa=True, higher_taxa=higher_taxa_name)
                                    sp_taxon_dict[species] = taxonomy
                                    #print(f"Corrected Species' ({species}) Higher Taxa Name Used to Find Taxa!\n")
                                except:
                                    sp_taxon_dict[species] = taxonomy
                                    #print(f'Finding Corrected Species Name for {species} Failed Due to Error Using Higher Taxa as Substitute.\n')  
                            else:
                                pass
                                #print(f'Corrected Species Name for {species} Used to Find Taxa!\n')
                        else:
                            pass
                            #print(f'Finding Corrected Species Name for {species} Failed.\n') 
                    except:
                        print(f'Finding Corrected Species Name for {species} Failed Due to Error.\n')  
            except:
                with open(taxon_file, 'w') as f:
                    json.dump(sp_taxon_dict, f, indent=4)  # indent for pretty formatting
                raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')

    return sp_taxon_dict

def ncbi_query_to_df(query_list, species_list, species_taxon_dict, email):
    
    """Converts a list of NCBI query results into a Pandas DataFrame.

    This function takes a list of Biopython SeqRecord objects (the results of NCBI queries) 
    and extracts relevant information such as accession numbers, taxonomic information, 
    gene descriptions, and protein sequences. It then organizes this data into a Pandas DataFrame.

    Args:
        query_list (list): A list of lists, where each inner list contains Biopython SeqRecord objects 
                           returned from NCBI queries for a specific species.
        species_list (list): A list of species names (strings) that were used in the NCBI queries.
        species_taxon_dict (dict): A dictionary containing taxonomic information for each species in 
                                   `species_list`, including synonyms. The keys are species names, 
                                   and the values are dictionaries with keys like "Phylum", "Subphylum", 
                                   "Class", and "Synonyms".
        email (str): Your email address, used for querying NCBI in the case of missing taxonomy.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted information from the NCBI query results. 
                          The DataFrame has the following columns:
                            - Accession: The NCBI accession number (version).
                            - Phylum: The phylum of the species.
                            - Subphylum: The subphylum of the species.
                            - Class: The class of the species.
                            - Genus: The genus of the species.
                            - Species: The species name (without the genus).
                            - Full_Species: The full species name (genus + species).
                            - Protein: The amino acid sequence of the protein.
                            - Gene_Description: A description of the gene.
                            - Species_Synonym_Used: The synonymous species name used in the query if a synonym was used, 
                                                    otherwise 'NA'.

    Raises:
        Exception: If a species query fails during the taxonomic lookup for an unknown species. This exception
                   is raised after 50 retries. This is highly unlikely, and likely indicates a back-end problem.
                   If you encounter this exception please save any work and restart the script.

    Notes:
      - The function handles cases where the species name in the NCBI record might be a synonym 
        of the name used in the original query.
      - It attempts to look up the taxonomic classification (Phylum, Subphylum, Class) 
        in the provided `species_taxon_dict`.
      - If a species is not found in `species_taxon_dict` and is not a synonym of any included species, it will attempt to fetch the taxonomic information from NCBI directly using the `get_species_taxonomy` function and the provided email. It will retry this process up to 50 times for each species if the initial attempt fails.
      - Duplicate entries (based on species name and protein sequence) are removed, keeping only the first occurrence.
    """
    # create empty lists
    Accession = []
    dna = []
    Phylum = []
    Order = []
    Class = []
    Genus = []
    Species = []
    gene_des = []
    version = []
    Protein = []
    full_sp_names = []
    Sp_syn_used = []
    
    # loop through the result list obtained from the NCBI search
    for query, sp in zip(query_list, species_list):
        # Get genus and speceis name seperately
        for seq in query:
            g_s_name = sp.split(' ', 1)
            # Search the dictionary of synonymous species names to see 
            # if this is the primary name or synonym.
            entry_spe_name = seq.annotations["organism"]
            
            l1 = len(entry_spe_name)
            l2 = len(sp)
            
            temp = entry_spe_name.split(' ')
            temp2 = sp.split(' ')
            if ((entry_spe_name == sp) or (entry_spe_name[:l1-1] == sp) or (entry_spe_name == sp[:l2-1])):
                found_w_synonym = False
                Phylum.append(species_taxon_dict[sp]["Phylum"])
                Class.append(species_taxon_dict[sp]["Class"])
                Order.append(species_taxon_dict[sp]["Order"])

            elif len(temp) == 3 or len(temp2) == 3:
                if (len(temp) == 3) and ((sp == str(temp[0]+' '+temp[1])) or (sp == str(temp[0]+' '+temp[2]))):
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                    Order.append(species_taxon_dict[sp]["Order"])
                elif (len(temp2) == 3) and ((entry_spe_name == str(temp2[0]+' '+temp2[1])) or (entry_spe_name == str(temp2[0]+' '+temp2[2]))):
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                    Order.append(species_taxon_dict[sp]["Order"])

                else:
                    g_s_name = entry_spe_name.split(' ', 1)
                    found_w_synonym = False
                    # If the species for this entry is not in the species list or synonyms dict then we will fetch the phylogeny now
                    queried = False
                    for x in range(50):
                        if queried == False:
                            try:        
                                temp_taxon = get_species_taxonomy(entry_spe_name, email)
                                queried = True
                            except:
                                pass
                        else:
                            pass
                    if queried == False:
                        raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')
                    Phylum.append(temp_taxon["Phylum"])
                    Class.append(temp_taxon["Class"])
                    Order.append(temp_taxon["Order"])
            else:
                if entry_spe_name in species_taxon_dict[sp]["Synonyms"]:
                    found_w_synonym = True
                    Phylum.append(species_taxon_dict[sp]["Phylum"])
                    Class.append(species_taxon_dict[sp]["Class"])
                    Order.append(species_taxon_dict[sp]["Order"])

                else:
                    g_s_name = entry_spe_name.split(' ', 1)
                    found_w_synonym = False
                    # If the species for this entry is not in the species list or synonyms dict then we will fetch the phylogeny now
                    queried = False
                    for x in range(50):
                        if queried == False:
                            try:        
                                temp_taxon = get_species_taxonomy(entry_spe_name, email)
                                queried = True
                            except:
                                pass
                        else:
                            pass
                    if queried == False:
                        raise Exception('Species query failed.\nThis is likely a back-end issue with the querying process\nSpecies Taxon Dictionary Saved. Please restart...')
                    Phylum.append(temp_taxon["Phylum"])
                    Class.append(temp_taxon["Class"])
                    Order.append(temp_taxon["Order"])

            if found_w_synonym == True:
                Sp_syn_used.append(entry_spe_name)
            else:
                Sp_syn_used.append('NA')

            if seq.seq:
                try:
                    #print(f'This is the DNA sequence: {seq.seq}')
                    dna_seq = seq.seq
                except:
                    dna_seq = ""
            else:
                dna_seq = ""
                
            pro_seq = ""
            # get and append protein sequence
            if seq.features:
                for feature in seq.features:
                    if feature.type == "CDS":
                        if "translation" in feature.qualifiers.keys():
                            pro_seq = feature.qualifiers['translation'][0]
                        
            # Append all meta data to corresponding lists
            Accession.append(str(seq.name))
            Genus.append(str(g_s_name[0]))
            Species.append(str(g_s_name[1]))
            full_sp_names.append(str(g_s_name[0]) + ' ' + str(g_s_name[1]))
            gene_des.append(str(seq.description))
            version.append(str(seq.id))
            dna.append(str(dna_seq))
            Protein.append(str(pro_seq))
            
    # create a dataframe for the information
    ncbi_q_df = pd.DataFrame(
        {'Accession': version,
        'Phylum': Phylum,
        'Class': Class,
        'Order': Order,
        'Genus': Genus,
        'Species': Species,
        'Full_Species': full_sp_names,
        'DNA': dna,
        'Protein': Protein,
        'Gene_Description': gene_des,
        'Species_Synonym_Used': Sp_syn_used
        })

    # Drop duplicates where the species names and protein sequences are the same...
    ncbi_q_df.drop_duplicates(subset=['Full_Species', 'Protein'],  keep='first', inplace=True)
    ncbi_q_df = ncbi_q_df.reset_index(drop=True)
    return ncbi_q_df
    


def ncbi_mine_seq_data(email, job_label='unnamed', out='unnamed', species_list=None, query=None, taxa_dictionary=None):
    
    """Mines NCBI for sequence data based on a list of species and a query.

    This function performs the following steps:
    1. Creates a directory to store results, labeled with the job name and current date/time.
    2. Saves the list of queried species to a text file.
    3. Constructs or loads a taxonomic dictionary for the species, including synonyms.
    4. Queries NCBI's databases for sequences matching the query for each species, considering synonyms.
    5. Formats the query results into a Pandas DataFrame.
    6. Saves the DataFrame to a CSV file.
    7. Saves the retrieved sequences to a FASTA file.
    8. Cleans the DataFrame to include only entries matching the input species list, saving cleaned and potential hits dataframes.
    9. Saves a cleaned FASTA file with sequences only from the input species list.
    10. Returns the resulting DataFrame (either cleaned or raw) and the report directory path.

    Args:
        email (str): Your email address, used for NCBI Entrez queries.
        job_label (str, optional): A label for the job. Defaults to 'unnamed'.
        out (str, optional): Output file name prefix. Defaults to 'unnamed'. If 'unnamed' it defaults to the job label.
        species_list (list, optional): A list of species names to query.
        query (str, optional): The search query for NCBI, e.g., "opsin AND Rhodopsin".
        taxa_dictionary (dict or str, optional): A dictionary or path to a JSON file containing a pre-built taxonomic dictionary. 
                                                 If None, a new dictionary will be constructed.

    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: The query results as a DataFrame (cleaned or original).
            - str: The path to the directory where results are saved.
    """
    
    print('Creating Job Directory\n')
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    report_dir = f'mnm_data/mnm_on_{job_label}_{dt_label}'
    os.makedirs(report_dir)
    
    print('Saving Species Query List to Text\n')
    with open(f"{report_dir}/species_queried.txt", "w") as f:
        for sp in species_list:
            f.write(str(sp) + "\n")
    
    if taxa_dictionary == None:
        # Create a taxonomic dictionary, including species synonyms, for all the species in the unique species list  
        print('Constructing Taxon Dictionary, Including Species Synonyms\n')
        # Check to see if an existing taxonomy dictionary already exists to save time.
        taxon_file = './data_sources/taxonomy/ncbi_taxon_dict.json'
        if os.path.isfile(taxon_file):
            try:
                with open(taxon_file, 'r') as f:
                    existing_taxon_dict = json.load(f)
            except FileNotFoundError:
                print(f"Error: File '{taxon_file}' not found or can't be loaded\n")
            print('Existing Taxon Dictionary Found! One Moment While We Update It...\n')
            species_taxon_dict = get_sp_taxon_dict(species_list = species_list, email = email, taxon_file = taxon_file, sp_taxon_dict = copy.deepcopy(existing_taxon_dict))
        else:
            existing_taxon_dict = {}
            species_taxon_dict = get_sp_taxon_dict(species_list = species_list, email = email, taxon_file = taxon_file)
        
        # Save the taxon dictionary if it doesn't yet exist or if it has been updated since being loaded 
        if (list(species_taxon_dict.keys()) >= list(existing_taxon_dict.keys())):         
            #print('Saving Updated Dictionary') 
            try:
                with open(taxon_file, 'w') as f:
                    json.dump(species_taxon_dict, f, indent=4)  # indent for pretty formatting
            except FileNotFoundError:
                print(f"Error: File '{taxon_file}' can't be saved...\n")

        print('Taxon Dictionary Complete!\n')
            
    else:
        if isinstance(taxa_dictionary, dict):
            species_taxon_dict = taxa_dictionary
        else:
            try:
                with open(taxa_dictionary, 'r') as f:
                    species_taxon_dict = json.load(f)
            except FileNotFoundError:
                raise Exception(f"Error: Taxon Dictionary File '{taxon_file}' not found or can't be loaded\n")
                
    # List to append query responses to
    query_list = []
    # make a progress bar for tracking query progression. Based on length of the species list
    print('Starting Queries to NCBI for DNA/Protein Sequences\n')
    i=0
    with progressbar.ProgressBar(max_value=len(species_list),style='BouncingBar') as bar:
        for species in species_list:
            
            try:
                temp = species.split(' ')
            except:
                raise Exception(f'Species Name Causing Error: {species}')
              
            if len(species_taxon_dict[species]['Synonyms']) > 0:
                sp_for_query = f'("{species}"[Organism] OR "{species}"[Title]'
                if (len(temp) == 3) and ('(' not in species) and ('.' not in species):
                    sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
                    sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism]  OR "{temp[0]} {temp[2]}"[Title]'
                    for syn in species_taxon_dict[species]['Synonyms']:
                        if (syn != species) and (f"{temp[0]} {temp[1]}" != syn) and (f"{temp[0]} {temp[2]}" != syn):
                            sp_for_query+= f' OR "{syn}"[Organism] OR {syn}"[Title]'
                            sp_for_query+=')'
                else:
                    for syn in species_taxon_dict[species]['Synonyms']:
                        if (syn != species):
                            sp_for_query+= f' OR "{syn}"[Organism] OR {syn}"[Title]'
                            sp_for_query+=')'
                            
            elif (len(temp) == 3) and ('(' not in species) and ('.' not in species):
                sp_for_query = f'("{species}"[Organism]'
                sp_for_query+= f' OR "{temp[0]} {temp[1]}"[Organism] OR "{temp[0]} {temp[1]}"[Title]'
                sp_for_query+= f' OR "{temp[0]} {temp[2]}"[Organism] OR "{temp[0]} {temp[2]}"[Title]'
                sp_for_query+=')'
                
            else:
                sp_for_query = f'"{species}"[Organism]'
            #print(f"{sp_for_query} AND {query}")
            queried = False
            for x in range(50):
                if queried == False:
                    try:            
                        NCBI_seq = ncbi_fetch(email=email, 
                                        term = f"{sp_for_query} AND {query}")
                        queried = True
                    except:
                        time.sleep(1)
                        pass
            if queried == False:
                print('Uh-oh, species query failed.\nThis is likely a back-end issue with the querying process\nIf this message continues to appear, please manually interrupt the qury process and restart...')
                    
            query_list.append(NCBI_seq)
            bar.update(i)
            i+=1
        bar.finish()
        
    print('NCBI Queries Complete!\nNow Extracting and Formatting Results For DataFrame...\n')
    ncbi_query_df = ncbi_query_to_df(query_list=query_list, species_list=species_list, species_taxon_dict=species_taxon_dict, email=email)
    
    if out == 'unnamed':
        out = job_label
    ncbi_query_df.to_csv(path_or_buf=f"./{report_dir}/{out.replace(' ','_')}_ncbi_q_data.csv", index=False)
    print('DataFrame Formatted and Saved to CSV file for future use :)\n')
    fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs.fasta'
    with open(fasta_file, 'w') as f:
        for id, seq in zip(ncbi_query_df['Accession'], ncbi_query_df['Protein']):
            f.write(f'>{id}\n{seq}\n')
    print('FASTA File Saved...\n')
    
    # Cleaning raw ncbi query df to keep only species entries that match our input sp list - returns the 'clean' df
    # Will also return a dataframe of the 'potential' hits - since it could be that the species that don't match is due to species mispelling or synonymous names
    ncbi_sp_hits = list(set(ncbi_query_df['Full_Species'].to_list()))
    #len(ncbi_sp_hits)
    intersection = list(set(ncbi_sp_hits) & set(species_list))
    #len(intersection)
    sp_no_hits  = list(set(species_list).symmetric_difference(intersection))
    #len(sp_no_hits)
    if len(sp_no_hits) > 0:
        print('Saving txt file with names of species that retrieved no results...\n')
        no_sp_hits_file = f'{report_dir}/species_w_no_hits.txt'
        with open(no_sp_hits_file, 'w') as f:
            for sp in sp_no_hits:
                f.write(f'{sp}\n')
    
    sp_rnd_hits  = list(set(ncbi_sp_hits).symmetric_difference(intersection))
    if len(sp_rnd_hits) > 0:
        print('Saving txt file with names of species that retrieved results but are NOT in submitted species list...\n')
        #len(sp_rnd_hits)
        rnd_sp_hits_file = f'{report_dir}/potential_species_hits.txt'
        with open(rnd_sp_hits_file, 'w') as f:
            for sp in sp_rnd_hits:
                f.write(f'{sp}\n')
        ncbi_query_df_cleaned = ncbi_query_df[~ncbi_query_df['Full_Species'].isin(sp_rnd_hits)]
        ncbi_query_df_potential_hits = ncbi_query_df[ncbi_query_df['Full_Species'].isin(sp_rnd_hits)]
        #ncbi_query_df_cleaned.shape
        #ncbi_query_df_potential_hits.shape
        print('Saving and returning cleaned dataframe with only species entries from species list...\n')
        ncbi_query_df_cleaned.to_csv(path_or_buf=f'{report_dir}/mnm_on_all_dbs_ncbi_q_data_cleaned.csv', index=False)
        print('Saving another dataframe with species that retrieved results but are NOT in submitted species list for further examination...\n')
        ncbi_query_df_potential_hits.to_csv(path_or_buf=f'{report_dir}/mnm_on_all_dbs_ncbi_q_potential_hits.csv', index=False)
        
        fasta_file = f'{report_dir}/mined_{out.replace(" ","_")}_seqs_cleaned.fasta'
        with open(fasta_file, 'w') as f:
            for id, seq in zip(ncbi_query_df_cleaned['Accession'], ncbi_query_df_cleaned['Protein']):
                f.write(f'>{id}\n{seq}\n')
        print('Clean FASTA File Saved...\n')
        return(ncbi_query_df_cleaned, report_dir)
        
    return(ncbi_query_df, report_dir)


from Bio import Entrez
import time
import random
import json

def ncbi_fetch_species(email, out='species_data', taxa=None, rank=None, limit=500, verbose=False):
    """
    Fetches species names from NCBI for a higher-order taxonomic group and 
    records taxonomic information for each species.

    Args:
        email (str): Your email address (required by NCBI).
        out (str): Output file name prefix (JSON file will be saved as out.json).
        taxa (str): The name of the higher-order taxonomic group (e.g., "Cestoda").
        rank (str): The taxonomic rank of the input taxa (e.g., "phylum", "class", "order").
        limit (int): The maximum number of species to fetch.

    Returns:
        dict: A dictionary where keys are species names and values are dictionaries 
              containing taxonomic information (Phylum, Class, Order, Family, etc.).
    """
    api_key = "1efb120056e1cea873ba8d85d6692abd5d09"
    Entrez.email = email
    ranks = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    rank = rank.capitalize()  # Standardize input for case-insensitivity
    primary_taxa_dict = {}
    
    try:
        # 1. Get the TaxID of the higher-order group
        print(f"Searching for TaxID for {taxa}...")
        search_handle = Entrez.esearch(db="taxonomy", term=f"{taxa}[{rank}]", retmax=1, api_key = api_key)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        if not search_record["IdList"]:
            print(f"Error: Could not find TaxID for {taxa} at rank {rank}.")
            return {}

        tax_id = search_record["IdList"][0]
        print(f"Found TaxID: {tax_id}")

        # 2. Find immediate children
        print(f"Finding direct children of {taxa}...")
        current_index = ranks.index(rank)
        child_rank = ranks[current_index + 1]
        child_search_handle = Entrez.esearch(db="taxonomy", term=f"(txid{tax_id}[Subtree] AND {child_rank}[Rank]) NOT species[Rank]", retmax=1000, api_key = api_key)
        child_search_record = Entrez.read(child_search_handle)
        child_search_handle.close()
        child_tax_ids = child_search_record["IdList"]

        if not child_tax_ids:
            print(f"No direct children found for {taxa}. Fetching species directly.")
            child_tax_ids = [tax_id]
        else:
            print(f"Found {len(child_tax_ids)} direct children.")
            
            child_fetch_handle = Entrez.efetch(db="taxonomy", id=child_tax_ids, retmode="xml", api_key = api_key)
            child_fetch_records = Entrez.read(child_fetch_handle)
            child_fetch_handle.close()
            for child_record in child_fetch_records:
                if verbose == True:
                    print(f"  - Child Taxon: {child_record['ScientificName']} (TaxID: {child_record['TaxId']})")

        # 3. Distribute the limit among the children
        species_per_child = limit // len(child_tax_ids)
        remainder = limit % len(child_tax_ids)

        # 4. Fetch species and their taxonomic info for each child
        for i, child_id in enumerate(child_tax_ids):
            num_to_fetch = species_per_child + (1 if i < remainder else 0)
            if verbose == True:
                print(f"\nFetching up to {num_to_fetch} species for child TaxID {child_id}...")

            queried = False
            for x in range(10):
                try:
                    if queried == False:
                        # 4.1 Get all species IDs under this child
                        species_search_handle = Entrez.esearch(db="taxonomy", term=f"txid{child_id}[Subtree] AND species[Rank]", retmax=100000, api_key = api_key)
                        species_search_record = Entrez.read(species_search_handle)
                        species_search_handle.close()
                        all_species_ids = species_search_record["IdList"]
                        queried = True
                except:
                    pass
                
            # 4.2 Sample species IDs
            if len(all_species_ids) <= num_to_fetch:
                sampled_species_ids = all_species_ids
            else:
                sampled_species_ids = random.sample(all_species_ids, num_to_fetch)

            # 4.3 Fetch species taxonomic information
            if len(sampled_species_ids) == 0:
                if verbose == True:
                    print("No species found.")
                continue
            
            queried = False
            for x in range(10):
                try:
                    if queried == False:            
                        species_fetch_handle = Entrez.efetch(db="taxonomy", id=sampled_species_ids, retmode="xml", api_key = api_key)
                        species_fetch_records = Entrez.read(species_fetch_handle)
                        #print(species_fetch_records)
                        species_fetch_handle.close()
                        queried = True
                except:
                    pass
                
            for species_record in species_fetch_records:
                species_name = species_record.get("ScientificName", species_record["TaxId"]) # Use TaxID if the name is not available
                if verbose == True:
                    print(f"  Processing: {species_name}")
                synonyms = []
                lineage = {}
                if "LineageEx" in species_record:
                    for taxon in species_record["LineageEx"]:
                        lineage[taxon["Rank"].capitalize()] = taxon["ScientificName"]
                    if ("OtherNames" in species_record):    
                        for name in species_record["OtherNames"]["Synonym"]:
                            if name not in synonyms:
                                if '(' in name:
                                    if ',' in name:
                                        rename = name.split('(')[0].strip()
                                        if (rename not in synonyms) and (rename != species_name): 
                                            synonyms.append(rename)
                                    else:
                                        rename = name.split('(')[1].strip()
                                        rename = rename.split(')')[0].strip() + ' ' + rename.split(')')[1].strip()
                                        if (rename not in synonyms) and (rename != species_name):
                                            synonyms.append(rename)
                                        rename2 = name.split('(')[0].strip() + rename.split(' ')[1]
                                        if (rename2 not in synonyms) and (rename2 != species_name):
                                            synonyms.append(rename2)
                                else:
                                    if ',' in name:
                                        rename = name.split(',')[0]
                                        rename = rename.split(' ')[0] + ' ' + rename.split(' ')[1]
                                        if (rename not in synonyms) and (rename != species_name):
                                            synonyms.append(rename)     
                                    else:
                                        synonyms.append(name)
                            else:
                                pass
                            
                primary_taxa_dict[species_name] = {
                    "Phylum": lineage.get("Phylum", "unknown"),
                    "Class": lineage.get("Class", "unknown"),
                    "Order": lineage.get("Order", "unknown"),
                    "Family": lineage.get("Family", "unknown"),
                    "Genus": lineage.get("Genus", "unknown"),
                    "Synonyms": synonyms,
                    "TaxId": species_record["TaxId"]
                }
            if verbose == True:
                print(f"Fetched and processed {len(species_fetch_records)} species.")
            time.sleep(1)
        # 5. Save the dictionary to a JSON file
        output_filename = f"{out}.json"
        with open(output_filename, "w") as outfile:
            json.dump(primary_taxa_dict, outfile, indent=4)

        print(f"\nTaxonomic data saved to {output_filename}")
        return primary_taxa_dict

    except Exception as e:
        print(f"An error occurred: {e}")
        return {}


from Bio import Entrez
import time
import random
import json

def ncbi_fetch_species(email, report_dir='taxonomy_data', out='species_data', taxa=None, rank=None, limit=500, verbose=False):
    """
    Fetches species names from NCBI for a higher-order taxonomic group and 
    records taxonomic information for each species.

    Args:
        email (str): Your email address (required by NCBI).
        report_dir (str): Report directory (folder) for storing taxonmy dictionary json files.
        out (str): Output file name prefix (JSON file will be saved as out.json).
        taxa (str): The name of the higher-order taxonomic group (e.g., "Cestoda").
        rank (str): The taxonomic rank of the input taxa (e.g., "phylum", "class", "order").
        limit (int): The maximum number of species to fetch.

    Returns:
        dict: A dictionary where keys are species names and values are dictionaries 
              containing taxonomic information (Phylum, Class, Order, Family, etc.).
    """
    api_key = "1efb120056e1cea873ba8d85d6692abd5d09"
    Entrez.email = email
    ranks = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    rank = rank.capitalize()  # Standardize input for case-insensitivity
    primary_taxa_dict = {}
    
    try:
        # 1. Get the TaxID of the higher-order group
        print(f"Searching for TaxID for {taxa}...")
        search_handle = Entrez.esearch(db="taxonomy", term=f"{taxa}[{rank}]", retmax=1, api_key = api_key)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        if not search_record["IdList"]:
            print(f"Error: Could not find TaxID for {taxa} at rank {rank}.")
            return {}

        tax_id = search_record["IdList"][0]
        print(f"Found TaxID: {tax_id}")

        # 2. Find immediate children
        print(f"Finding direct children of {taxa}...")
        current_index = ranks.index(rank)
        child_rank = ranks[current_index + 1]
        child_search_handle = Entrez.esearch(db="taxonomy", term=f"(txid{tax_id}[Subtree] AND {child_rank}[Rank]) NOT species[Rank]", retmax=1000, api_key = api_key)
        child_search_record = Entrez.read(child_search_handle)
        child_search_handle.close()
        child_tax_ids = child_search_record["IdList"]

        if not child_tax_ids:
            print(f"No direct children found for {taxa}. Fetching species directly.")
            child_tax_ids = [tax_id]
        else:
            print(f"Found {len(child_tax_ids)} direct children.")
            
            child_fetch_handle = Entrez.efetch(db="taxonomy", id=child_tax_ids, retmode="xml", api_key = api_key)
            child_fetch_records = Entrez.read(child_fetch_handle)
            child_fetch_handle.close()
            for child_record in child_fetch_records:
                print(f"  - Child Taxon: {child_record['ScientificName']} (TaxID: {child_record['TaxId']})")

        # 3. Distribute the limit among the children
        species_per_child = limit // len(child_tax_ids)
        remainder = limit % len(child_tax_ids)

        # 4. Fetch species and their taxonomic info for each child
        for i, child_id in enumerate(child_tax_ids):
            num_to_fetch = species_per_child + (1 if i < remainder else 0)
            print(f"\nFetching up to {num_to_fetch} species for child TaxID {child_id}...")

            queried = False
            for x in range(10):
                try:
                    if queried == False:
                        # 4.1 Get all species IDs under this child
                        species_search_handle = Entrez.esearch(db="taxonomy", term=f"txid{child_id}[Subtree] AND species[Rank]", retmax=100000, api_key = api_key)
                        species_search_record = Entrez.read(species_search_handle)
                        species_search_handle.close()
                        all_species_ids = species_search_record["IdList"]
                        queried = True
                except:
                    pass
                
            # 4.2 Sample species IDs
            if len(all_species_ids) <= num_to_fetch:
                sampled_species_ids = all_species_ids
            else:
                sampled_species_ids = random.sample(all_species_ids, num_to_fetch, )

            # 4.3 Fetch species taxonomic information
            if len(sampled_species_ids) == 0:
                print("No species found.")
                continue
            
            queried = False
            for x in range(10):
                try:
                    if queried == False:            
                        species_fetch_handle = Entrez.efetch(db="taxonomy", id=sampled_species_ids, retmode="xml", api_key = api_key)
                        species_fetch_records = Entrez.read(species_fetch_handle)
                        print(species_fetch_records)
                        species_fetch_handle.close()
                        queried = True
                except:
                    pass
                
            for species_record in species_fetch_records:
                species_name = species_record.get("ScientificName", species_record["TaxId"]) # Use TaxID if the name is not available
                print(f"  Processing: {species_name}")
                synonyms = []
                lineage = {}
                if "LineageEx" in species_record:
                    for taxon in species_record["LineageEx"]:
                        lineage[taxon["Rank"].capitalize()] = taxon["ScientificName"]
                    if ("OtherNames" in species_record):    
                        for name in species_record["OtherNames"]["Synonym"]:
                            if name not in synonyms:
                                if '(' in name:
                                    if ',' in name:
                                        rename = name.split('(')[0].strip()
                                        if (rename not in synonyms) and (rename != species_name): 
                                            synonyms.append(rename)
                                    else:
                                        rename = name.split('(')[1].strip()
                                        rename = rename.split(')')[0].strip() + ' ' + rename.split(')')[1].strip()
                                        if (rename not in synonyms) and (rename != species_name):
                                            synonyms.append(rename)
                                        rename2 = name.split('(')[0].strip() + rename.split(' ')[1]
                                        if (rename2 not in synonyms) and (rename2 != species_name):
                                            synonyms.append(rename2)
                                else:
                                    if ',' in name:
                                        rename = name.split(',')[0]
                                        rename = rename.split(' ')[0] + ' ' + rename.split(' ')[1]
                                        if (rename not in synonyms) and (rename != species_name):
                                            synonyms.append(rename)     
                                    else:
                                        synonyms.append(name)
                            else:
                                pass
                            
                primary_taxa_dict[species_name] = {
                    "Phylum": lineage.get("Phylum", "unknown"),
                    "Class": lineage.get("Class", "unknown"),
                    "Order": lineage.get("Order", "unknown"),
                    "Family": lineage.get("Family", "unknown"),
                    "Genus": lineage.get("Genus", "unknown"),
                    "Synonyms": synonyms,
                    "TaxId": species_record["TaxId"]
                }

            print(f"Fetched and processed {len(species_fetch_records)} species.")
            time.sleep(1)
        # 5. Save the dictionary to a JSON file
        output_filename = f"./{report_dir}/{out}.json"
        with open(output_filename, "w") as outfile:
            json.dump(primary_taxa_dict, outfile, indent=4)

        print(f"\nTaxonomic data saved to {output_filename}")
        return primary_taxa_dict

    except Exception as e:
        print(f"An error occurred: {e}")
        return {}
