import subprocess  # For running external programs like tblastx
import datetime    # For timestamping output filenames
import os          # For interacting with the filesystem

# ==========================
# RUN BLASTX
# ==========================
def run_blast_query(blasttyp, QUERY_FILE, DB_FILE, OUTPUT_FILE, EVALUE, OUTFMT):
    # Construct blastx command as a list of arguments
    blast_cmd = [
        blasttyp,             # BLAST tool to run
        "-query", QUERY_FILE,  # Your transcriptome FASTA file
        "-db", DB_FILE,        # Database to search against
        "-out", OUTPUT_FILE,   # Where to save the results
        "-evalue", EVALUE,     # E-value threshold
        "-outfmt", OUTFMT      # Output format (e.g., tabular format)
    ]

    print(f"Running {blasttyp}...")
    result = subprocess.run(blast_cmd, capture_output=True, text=True)

    if result.stderr:
        print(f"Error running {blasttyp}:", result.stderr)
    else:
        print(f"{blasttyp} search completed successfully!")
        print(f"Results saved in: {OUTPUT_FILE}")

    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, "r") as f:
            lines = f.readlines()
            print("\n=== Preview of Results ===")
            print("".join(lines[:10]))
    else:
        print("Error: No output file found. Check if blastx ran correctly.")

# ==========================
# FORMAT BLAST DATABASE (protein)
# ==========================
def format_blast_db(db_file):
    print(f"Formatting BLAST DB for: {db_file}")
    db_base = os.path.splitext(db_file)[0]  # Strip .fasta
    result = subprocess.run([
        "makeblastdb",
        "-in", db_file,
        "-dbtype", "prot",
        "-out", db_base  # Explicitly set base name of DB
    ], capture_output=True, text=True)

    if result.stderr:
        print("makeblastdb error:", result.stderr)
    else:
        print("Database formatted successfully.")

# ==========================
# MAIN
# ==========================
def main():
    print("Script started!")
    # User settings
    DB_folder = '/home/local/ADS/nicoleharris/labdata/users/McKinley/Nicole_Harris/blast_dbs/opsin_dbs'
    QUERY_folder = "/home/local/ADS/nicoleharris/labdata/ostracod_seqData/aa"
    job_name = 'transx'
    results_folder = '/home/local/ADS/nicoleharris/labdata/users/McKinley/Nicole_Harris/opsin_blast_results'

    if not os.path.isdir(results_folder):
        os.makedirs(results_folder, exist_ok=True)

    # blastx parameters
    evalue = "1e-5"
    outfmt = "6"
    blasttyp = 'blastp'

    print("Checking DB folder:", DB_folder)
    print("Checking QUERY folder:", QUERY_folder)

    for db in os.listdir(DB_folder):
        if db.endswith('.fasta'):
            print("Found DB file:", db)
            for query in os.listdir(QUERY_folder):
                if query.lower().endswith(('.fasta', '.fa', '.faa', '.pep')):
                    print("Found QUERY file:", query)

                    query_path = os.path.join(QUERY_folder, query)
                    db_path = os.path.join(DB_folder, db)

                    # Format DB and get base path
                    format_blast_db(db_path)
                    db_base = os.path.splitext(db_path)[0]

                    # Output file
                    query_name = os.path.splitext(query)[0]
                    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                    output_file = os.path.join(results_folder, f"{job_name}_{blasttyp}_{query_name}_{time_stamp}.txt")

                    run_blast_query(blasttyp, query_path, db_base, output_file, evalue, outfmt)

if __name__ == '__main__':
    main()
    