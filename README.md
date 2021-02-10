# HbxFinder

HbxFinder is a pipeline for identification and characterisation of homeobox genes in unannotated genomes. Using raw data from [HomeoDB](http://homeodb.zoo.ox.ac.uk/) by [Zhong & Holland 2011](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1525-142X.2011.00513.x), this tool uses BLAST and [MMseqs2](https://github.com/soedinglab/MMseqs2) similarity searches to annotate gene families containing homeodomains. Currently, this tool is written for use on high quality arthropod genomes.

## Requirements

* python3
* BLAST+ suite - [install for command line](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* MMseqs - [install](https://github.com/soedinglab/MMseqs2#installation)

## Usage

HbxFinder is a collection of python script that is currently run with a string of commands.

### Step 1: Run initial broad tBLASTn search

`python blast_run.py --path </path/to/genome/assemblies/> `

- Required input:
  - **`directory`** First change dir with `cd HbxFinder/hbxfinder/blast_search`. This puts you in the correct directory to run the blast script.
  - **`genomes`** Path to raw genome assemblies you wish to search.
- Output: Blast output files in output format 6. The most important file for the next step is present at `recip_blast/genome_recipBlast.fa`

Following the initial tBLASTn search, a reciprocal BLASTx search is carried out by running:

`python recip_blast.py`

Followed by:

`python summarise_Hbx.py --taxa <reciprocal blastoutput file> --gene <specify class of hbx gene>`

This produces an output file for each homeobox gene whic is parsed in Step 2.

### Step 2: Run	second more sensitive sequence similarity search with MMseqs

The second similarity search uses MMseqs to carry out a more sensitive search of the regions of the genome containing homeobox genes. You must specify a particular class of homeobox gene you wish to annotated e.g. HOX.

`python data_blast.py --gene <specify class of hbx gene> --path </path/to/genome/assemblies/> `


### Step 3: Get sequence data for homeobox genes

`python get_seq_data.py --gene <specify class of hbx gene> `