import subprocess 
import datetime    
import os          
import pandas as pd 

# Define standard column names for BLAST output format 6
BLAST_COLUMN_NAMES = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

# ==========================
# RUN BLASTP/BLASTX
# ==========================
def run_blast_query(blasttyp, QUERY_FILE, DB_BASE_PATH, OUTPUT_FILE, EVALUE, OUTFMT, column_names):
    """
    Runs a BLAST search and returns the results as a pandas DataFrame.
    DB_BASE_PATH should be the base path for the formatted BLAST database (without .fasta or extensions like .pin).
    """
    blast_cmd = [
        blasttyp,
        "-query", QUERY_FILE,
        "-db", DB_BASE_PATH,    # Database base path to search against
        "-out", OUTPUT_FILE,    # Where to save the raw results
        "-evalue", EVALUE,      # E-value threshold
        "-outfmt", OUTFMT       # Output format (e.g., tabular format 6)
    ]

    #print(f"Running {blasttyp} with query '{os.path.basename(QUERY_FILE)}' against DB '{os.path.basename(DB_BASE_PATH)}'...")
    result = subprocess.run(blast_cmd, capture_output=True, text=True)

    if result.stderr:
        print(f"  Error running {blasttyp}: {result.stderr.strip()}")
        return pd.DataFrame()  # Return empty DataFrame on BLAST error

    #print(f"  {blasttyp} search completed. Raw results intended for: {OUTPUT_FILE}")

    if os.path.exists(OUTPUT_FILE) and os.path.getsize(OUTPUT_FILE) > 0:
        try:
            blast_df = pd.read_csv(OUTPUT_FILE, sep='\t', header=None, names=column_names)
            #print(f"  Successfully read {len(blast_df)} BLAST hits into DataFrame from {OUTPUT_FILE}.")
            return blast_df
        except pd.errors.EmptyDataError:
            print(f"  Warning: BLAST output file {OUTPUT_FILE} is empty.")
            return pd.DataFrame()
        except Exception as e:
            print(f"  Error reading BLAST output file {OUTPUT_FILE} into DataFrame: {e}")
            return pd.DataFrame()
    else:
        if not os.path.exists(OUTPUT_FILE):
            print(f"  Error: Output file {OUTPUT_FILE} was not created by BLAST.")
        elif os.path.getsize(OUTPUT_FILE) == 0:
            print(f"  Warning: BLAST output file {OUTPUT_FILE} is empty (no hits found or error).")
        return pd.DataFrame()

# ==========================
# FORMAT BLAST DATABASE (protein)
# ==========================
def format_blast_db(db_file_path):
    """
    Formats a FASTA file into a BLAST database.
    """
    print(f"Formatting BLAST DB for: {db_file_path}")
    db_base = os.path.splitext(db_file_path)[0]  # Strip .fasta to get base name for DB files
    
    # Check if database files already exist to avoid re-formatting if not necessary
    # makeblastdb typically creates multiple files (e.g., .pin, .psq, .phr)
    # A simple check for one of them, e.g., .pin or .psq depending on version/type
    # For protein dbs, .psq is common.
    required_db_files_exist = os.path.exists(db_base + ".pdb") or os.path.exists(db_base + ".pin")
    
    if required_db_files_exist:
        print(f"  Database files for '{db_base}' seem to exist. Skipping formatting.")
        # You might add a check here to see if the source .fasta is newer than the DB files
        # and force re-formatting if it is. For simplicity, this is omitted here.
        return True # Indicate success or that formatting is not needed

    result = subprocess.run([
        "makeblastdb",
        "-in", db_file_path,
        "-dbtype", "prot",
        "-out", db_base  # Explicitly set base name of DB
    ], capture_output=True, text=True)

    if result.stderr:
        print(f"  makeblastdb error for {db_file_path}: {result.stderr.strip()}")
        return False
    else:
        print(f"  Database '{db_base}' formatted successfully.")
        return True

# ==========================
# MAIN
# ==========================
def main():
    print("Script started!")
    # User settings
    DB_folder = '/home/local/ADS/nicoleharris/labdata/users/McKinley/Nicole_Harris/blast_dbs/opsin_dbs'
    QUERY_folder = "/home/local/ADS/nicoleharris/labdata/ostracod_seqData/aa"
    job_name = 'transx' # This can represent your protein family or project name
    results_folder = '/home/local/ADS/nicoleharris/labdata/users/McKinley/Nicole_Harris/opsin_blast_results'

    if not os.path.isdir(results_folder):
        os.makedirs(results_folder, exist_ok=True)
        print(f"Created results folder: {results_folder}")

    # blastp parameters
    evalue = "1e-5"
    outfmt = "6" # Tabular format, essential for pandas import
    blasttyp = 'blastp' # Ensure this matches the type of query and DB

    print(f"DB folder: {DB_folder}")
    print(f"QUERY folder: {QUERY_folder}")
    print(f"Results will be saved in: {results_folder}")
    print(f"Job name / Protein family: {job_name}")
    print("---")

    for db_filename in os.listdir(DB_folder):
        if db_filename.endswith('.fasta') or db_filename.endswith('.faa'):
            print(f"\nProcessing Database File: {db_filename}")
            db_path_full = os.path.join(DB_folder, db_filename)
            db_name_base = os.path.splitext(db_filename)[0].replace('.', '_') # For naming output files

            # Format DB (once per database file)
            if not format_blast_db(db_path_full):
                print(f"  Skipping database {db_filename} due to formatting error.")
                continue # Move to the next database file
            
            # This is the base path that BLAST tools will use with the -db argument
            db_blast_base_path = os.path.splitext(db_path_full)[0]

            all_best_hits_for_this_db_list = [] # Accumulate best hits for THIS DB from all query files

            for query_filename in os.listdir(QUERY_folder):
                if query_filename.lower().endswith(('.fasta', '.fa', '.faa', '.pep')):
                    print(f"  Processing Query File: '{query_filename}' against DB: '{db_filename}'")
                    query_path_full = os.path.join(QUERY_folder, query_filename)
                    query_name_base = os.path.splitext(query_filename)[0]

                    # Temporary output file for this specific BLAST run's raw results
                    time_stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                    temp_raw_blast_output_file = os.path.join(
                        results_folder,
                        f"temp_blast_{job_name}_{query_name_base}_vs_{db_name_base}_{time_stamp}.tsv"
                    )

                    # Run BLAST and get DataFrame
                    blast_results_df = run_blast_query(
                        blasttyp, query_path_full, db_blast_base_path,
                        temp_raw_blast_output_file, evalue, outfmt, BLAST_COLUMN_NAMES
                    )

                    if blast_results_df is not None and not blast_results_df.empty:
                        # Ensure 'pident' is numeric for correct identification of maximum
                        blast_results_df['pident'] = pd.to_numeric(blast_results_df['pident'], errors='coerce')
                        blast_results_df.dropna(subset=['pident'], inplace=True) # Remove rows if pident wasn't numeric

                        if not blast_results_df.empty:
                            # Filter to keep only the hit with the highest 'pident' for each 'qseqid'
                            # Using .copy() to avoid SettingWithCopyWarning
                            idx = blast_results_df.groupby('qseqid', group_keys=False)['pident'].idxmax()
                            best_hits_df = blast_results_df.loc[idx].copy()
                            
                            best_hits_df['source_transcriptome_file'] = query_filename # Add source file info
                            all_best_hits_for_this_db_list.append(best_hits_df)
                            print(f"    Found {len(best_hits_df)} best hits from '{query_filename}'.")
                        else:
                            print(f"    No valid numeric 'pident' values in results for '{query_filename}' against '{db_filename}'.")
                        
                        # Optional: Remove the temporary raw BLAST output file if it's no longer needed
                        # if os.path.exists(temp_raw_blast_output_file):
                        #     os.remove(temp_raw_blast_output_file)
                        #     print(f"    Removed temporary file: {temp_raw_blast_output_file}")

                    else:
                        print(f"    No results or error processing BLAST output for '{query_filename}' against '{db_filename}'.")
                        # Optional: remove empty temp file if created
                        # if os.path.exists(temp_raw_blast_output_file) and os.path.getsize(temp_raw_blast_output_file) == 0:
                        #    os.remove(temp_raw_blast_output_file)


            # After processing all query files for the current database
            if all_best_hits_for_this_db_list:
                final_compiled_df_for_db = pd.concat(all_best_hits_for_this_db_list, ignore_index=True)
                
                # Final compiled CSV filename: jobName_dbName_bestHits.csv
                compiled_csv_filename = os.path.join(results_folder, f"{job_name}_{db_name_base}_best_hits_compiled.csv")
                final_compiled_df_for_db.to_csv(compiled_csv_filename, index=False)
                
                print(f"\n  Successfully compiled {len(final_compiled_df_for_db)} total best hits for database '{db_filename}'.")
                print(f"  Results saved to: {compiled_csv_filename}")
            else:
                print(f"\n  No best hits found for any query files against database '{db_filename}'.")
            print("---")

    print("\nScript finished!")

if __name__ == '__main__':
    main()
    