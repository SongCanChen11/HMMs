The custom script `hmmsearch.sh` aims to retrieve sulfur-cycling proteins from the GTDB genome database (https://gtdb.ecogenomic.org/) using hidden markov models (HMMs). The HMMs for sulfur-cycling proteins (n = 116) are compiled in the zipped file named `HMMs.tar.gz`. The description for each HMM is detailed in the `hmm_info.xlsx`.  A subset of sequences (~30 M unzipped; `GTDB_r95_demo.faa.gz`) from the full GTDB r95 database (~40 G) was used as a demo dataset.


# System requirement
## Hardware requirements
The `hmmsearch.sh` script requires a standard computer with >4G RAM and >5 CPUs to run the Demo dataset, and a server with >10G RAM and >20 CPUs to run the full GTDB r95 dataset.

## Software requirements
### OS Requirements
This package is supported for *macOS* and *Linux*. The script has been tested on the following systems:
+ macOS: Mojave (10.14.6)
+ Linux: Ubuntu 18.04.3 LTS

### Dependencies
The `hmmsearch.sh` mainly depends on the `hmmsearch` command in the HMMER package (v3.2; http://hmmer.org/).
The `parse_hmmscan.py` mainly depends on the pandas package of python.


# Installation Guide:

### install HMMER and BLAST+ using conda
The installation of hmmer/BLAST+ takes less than 5 minutes with a good connection to the conda mirror.
```
conda install hmmer
conda install blast
pip install pandas
```


# Run the Demo
The execution of the `hmmsearch.sh` script on the Demo dataset (GTDB_r95_demo.faa) takes ~30 min on a computer with 5 CPUs. The script will generate tabulated output for each sulfur-cycling gene in the Res directory. To run `hmmsearch.sh` on the full GTDB database, the `GTDB_r95_demo.faa` in the script should be replaced with the path to the dataset download from GTDB (https://gtdb.ecogenomic.org/).  
```
# unzip the HMMs.tar.gz file
tar -xzvf HMMs.tar.gz

# construct binary compressed datafiles for hmmscan
cat HMMs/*.curated >hmmdb.hmm
hmmpress hmmdb.hmm

# unzip the demo dataset
# create blastp database of demo dataset for sequence extraction
gunzip GTDB_r95_demo.faa.gz
makeblastdb -dbtype prot -out GTDB_r95_demo -in GTDB_r95_demo.faa -parse_seqids  


# make a result directory
mkdir Res

# Execute the hmmsearch of marker genes against the demo dataset
source hmmsearch.sh

# merge all sulfur-cycling homologs into a single result file named '01_merged_res.uniq.id' (n = 48537)
rm 01_merged_res.id
for file in `ls Res/*.ga.out`;
do 
	grep -v '^#' $file| awk '{print $1}' >>01_merged_res.id
done
sort 01_merged_res.id |sort|uniq >01_merged_res.uniq.id
wc 01_merged_res.uniq.id

# extract the target sequences from the demo dataset using blastcmd
blastdbcmd -db GTDB_r95_demo -dbtype prot -entry_batch 01_merged_res.uniq.id >01_merged_res.uniq.faa


# Annotate all candidates using hmmscan
hmmscan --cut_ga --tblout 01_merged_res.hmmscan.txt --cpu 20 hmmdb.hmm 01_merged_res.uniq.faa &>/dev/null &

# Parse the results of hmmscan; best hit was kept as annotation of the sequence
python parse_hmmscan.py -i 01_merged_res.hmmscan.txt -f hmm_info.xlsx -o 01_merged_res.anno.txt

```

# License
This project is covered under the **Apache 2.0 License**.


