Dataset containing available ostracod genomic/transcriptomic data put together in fall 2024 by Marie Drabkova


right now the dataset is located on osiris/labdata/ostracod_seqData
included is this readme file, and a spreadsheet with file information which is also available on a following link:
https://docs.google.com/spreadsheets/d/1Qr0cGsnSZRjhH4SzAJMwaAOmzI_39EvyoyMZwUllhxQ/edit?usp=sharing 


Most of the transcriptomes were already on Osiris or Pod servers and were processed by Depeng Li or other people in the lab. Some were originally downloaded from ncbi sra and have accession numbers (short reads archive), some came from Emily Elis’ paper (10.5281/zenodo.5498416) 
https://zenodo.org/records/5498417 
https://datadryad.org/stash/dataset/doi:10.6076/D1PC7K 
https://academic.oup.com/sysbio/article/72/2/264/6671865#408728285 and some from Todd’s older paper https://doi.org/10.5061/dryad.tb40v 
https://zenodo.org/records/5025796 


Some transcriptome or genome data came from matedb2, a recently published set of decently assembled and processed transcriptomes and genomes across Metazoa (also a paper). https://github.com/MetazoaPhylogenomicsLab/MATEdb2/ 


Some of the Depeng’s transcriptomes has not been published yet (in pink in the spreadsheet)


I also added two Jimmorinia transcriptomes assembled by Cory.


And two ostracods that were part of Rebecca Varney’s run of genome skimming (low coverage genome sequencing). These are called cave cod and hydroterm cod. I trimmed the reads and assembled them in platanus/megahit (more info on that in /labdata/genomeSkim01). Unfortunately these are not going to be much useful for genomic analyses, busco scores are really low, but they might be used for mitochondrial genome analyses.  


In the first tab of the spreadsheet you should be able to see the names of files and their respective sizes. On the second tab you can find more info about samples. I also added busco scores (evaluation of completeness of transcriptome, see https://busco.ezlab.org/, arthropoda set). Results are in busco_results folder (zipped files: busco_ostracods.tar.gz transcriptomesDL_busco.tar.gz, busco_transcriptomesEE.tar.gz). 


I usually use protein files. If you need nucleotide files for the corresponding protein files you might need to go back to the sources (like matedb2 and zenodo/dryad).


Processing of transcriptomes is pretty standard and straightforward, it is basically the same for all included data. Download reads (from database prefetch, fasterq-dump, or a sequencing center) > trimm (trimmomatic)> assemble (Trinity) > predict proteins (Transdecoder).
  

legend: dpl_np means from Depeng Li not published, RV Rebecca Varney, others are I hope clear. It is just a theoretical tree, not a real one!


I am also attaching busco script for pod as an example of pod script: SLURM_busco_tr.sh
I use largemem queue (option -p) on pod for this, it means it can use as much memory (RAM) as needed. Sometimes the queue is shorter than for normal batch jobs, sometimes you have to wait longer. 
I installed busco in environment with micromamba (anaconda variation), that needs to be activated first.
Plus it did not like direct download of the arthropod reference file, so I downloaded it before the run manually and ran this in the offline version.
Then I loop over files I need to process. 
Pod sometimes has an issue with the number of files on the system. To avoid that you can work on bigscratch (create a folder there and work there), but be careful about backup.


#!/bin/bash  -l
#SBATCH --nodes=1 --ntasks-per-node 20 -p largemem
#SBATCH --time=72:00:00


export MAMBA_ROOT_PREFIX=/home/mariedrabek/progs/micromamba/
eval "$(/home/mariedrabek/progs/micromamba/bin/micromamba shell hook -s posix)"
micromamba activate busco


cd /bigscratch/mariedrabek/busco


list=$(ls -1 /home/mariedrabek/matedbSingularity/transcriptomes_sra/|grep mod|grep pep)
for f in $list; do
 busco -i /home/mariedrabek/matedbSingularity/transcriptomes_sra/$f -m proteins -l arthropoda_odb10 -c 20 -o tr_busco_$f -f --offline
done