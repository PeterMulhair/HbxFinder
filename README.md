# HbxFinder

<div align="center">
<p align="center">
<img src="https://github.com/PeterMulhair/HbxFinder/blob/master/hbxfinder_logo.png" width="500" height="150">
</p>
</div>

HbxFinder is a pipeline for identification and characterisation of homeobox genes in unannotated genomes. Using raw data from [HomeoDB](http://homeodb.zoo.ox.ac.uk/) by [Zhong & Holland 2011](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1525-142X.2011.00513.x), this tool uses BLAST and [MMseqs2](https://github.com/soedinglab/MMseqs2) similarity searches to annotate gene families containing homeodomains. This tool is optimised to be used on either arthropod or vertebrate genomes.


## Requirements

* python3
* BLAST+ suite - [install for command line](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* MMseqs - [install](https://github.com/soedinglab/MMseqs2#installation)
* Emboss sixpack

## Usage

HbxFinder is a collection of python scripts that is currently run with a string of commands.

To run, download this repo locally using `git clone https://github.com/PeterMulhair/HbxFinder.git`

### Step 1: Run initial broad tBLASTn search

#### 1.1

`python blast_run.py --path </path/to/genome/assemblies/> --group <group of organisms to search>` (`--group` is given as either `invertebrate` or `vertebrate` depending on what group of animals you are searching)

- Required:
  - **`directory`** First change dir with `cd HbxFinder/hbxfinder/blast_search`. This puts you in the correct directory to run the blast script.
  - **`genomes`** Path to a folder containing the raw genome assemblies you wish to search. It is preferred if you provide the full path name to where genomes are stored.
- Output: Blast output files in output format 6. The most important file for the next step is present at `recip_blast/blast_parsed_output.fasta`

#### 1.2

Following the initial tBLASTn search, a reciprocal BLASTx search (to identify the hbx gene) is carried out followed by a parsing script (which orders the non-overlapping hbx genes) by running:

```
$ cd recip_blast/

$ python recip_blast.py --infile <output from blast script i.e. blast_parsed_output.fasta> --outfile <output file name> --group <group of organisms to search i.e. invertebrate or vertebrate>

$ python summarise_Hbx.py --taxa <reciprocal blastoutput file i.e. file ending in *_hbx.tsv> --gene <specify class of hbx gene> --group <invertebrate or vertebrate>

```

This produces an output file for each homeobox gene class (in the output dir `hbx_contigs`) which is parsed in Step 2.

### Step 2: Run	second more sensitive sequence similarity search with MMseqs

#### 2.1

The second similarity search uses MMseqs to carry out a more sensitive search of the regions of the genome containing homeobox genes. You must specify a particular class of homeobox gene you wish to annotated e.g. HOX.

`python data_blast.py --gene <specify class of hbx gene> --path </path/to/genome/assemblies/> --group <invertebrate or vertebrate>`


- Required:
  - **`directory`** First change dir with `cd HbxFinder/hbxfinder/second_search`. This puts you in the correct directory to run the blast script.
  - **`genomes`** Path to raw genome assemblies you wish to search. Again, it is preferred if you provide the full path name to where genomes are stored.
  - **`gene`** Name of homeobox gene class you wish to search for eg. HOX

#### 2.2

This is followed again by a reciprocal BLAST script along with a parsing script.

```
$ cd recip_blast/

$ python recip_blast.py --gene <specify class of hbx gene> --group <invertebrate or vertebrate>

$ python summarise_Hbx.py --taxa species_hbx/<species reciprocal blastoutput file> --gene <specify class of hbx gene> --path <path to genome file> --group <invertebrate or vertebrate>

```

The final output files are found in `hbx_clusters/`

### Step 3: Get sequence data for homeobox genes

The last step outputs fasta files of the nucleotide and amino acid sequences for the homeodomain of each homeobox gene class.

`python get_seq_data.py --gene <specify class of hbx gene> `


- Required:
  - **`directory`** First change dir with `cd HbxFinder/hbxfinder/second_search/HD_seq_out`
  - **`gene`** Name of homeobox gene class you wish to parse eg. HOX


---

## NOTE

- This pipeline only annotates the homeodomain sequence rather than the full open reading frame.
- It is recommended to manually inspect and edit the cluster files as well as the output fasta files to ensure accuracy in the order and content of the homeobox genes eg. if there is an intron in the homeodomain sequence, that gene will be annotated in two parts - you can edit this manually by checking the size of the homeodomain sequences.
- To get most accurate output, particularly for homeobox genes grouped in a cluster, it is recommended to run this script on high quality genome assemblies.
- If you have any questions or run into any issues please feel free to contact me by raising an issue or emailing me at `peter.mulhair@zoo.ox.ac.uk`