# Proteomics Simulator

This is simulator for proteomics prosidures for the detection of certain peptides and their mutated (edited, by RNA editing) versopns.
It consists of a few modules for the creation of a DB containing all the peptides and their edited versions that are possible in a proteomics procedure with certain parameters.

It also have a few modules for invoking [MaxQuant](https://www.maxquant.org/) with the created DB on multiple real proteomics samples in order to detect and quantify the expression of certain genes in their native and edited versions.

There is also a module that will help you create a control group of "randomly modified proteins" in order to asses the real editing signal in proteomics samples and quantify the SNR,


## Prerequisits

1. Python 3.6 or higher. including [Biopython](https://biopython.org/docs/1.75/api/Bio.html)

2. [MaxQuant](https://www.maxquant.org/), only mandatory for the analysis of proteoms in real proteomics samples, but not for the creation of the peptides data-based


## Data inputs for the Simulator

#### 1. Transcriptome:

A fasta file containing the coding regions of the mRNA of different genes
Each record Id in the fasta file will be considered as the transctipt name 

#### 2. Editome

A table containing coordinates (relative to the coding region + genomic coordinates) of different editing events

| transcript_name;gene_name | mismatch coordinate (base0) wrt coding transctip | mismatch coordinate (base1) wrt coding transctip | mismatch type | chromosome | mismatch coordinate (base0) wrt chromosome | mismatch coordinate (base1) wrt chromosome | strand direction of transctiprt | amino-acid swap | editing level (in case you examining RNA editing sites |
|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|
|NM_000047.2;ARSE        |1602    |1603    |TC      |chrX    |2934998 |2934999 |-       |YH      |0.190372|
|NM_000090.3;COL3A1      |4300    |4301    |GT      |chr2    |189011673       |189011674       |+       |RL      |0.0393992|
|NM_000113.2;TOR1A       |753     |754     |TC      |chr9    |129814216       |129814217       |-       |FL      |0.0163923|


You can also have consensus transcripts from some published genome and coordinates of mutation within that genome for which you would like to asses the protein expression.

Command for running the simulator.

`python proteomics_simulator_full_pipline.py -i_fasta <input_mRNAs> -i_sites <editing_sites_file> -mc 3 -min_aa 6 -max_m 5200`


This command uses many modules in a few stages. you can use them separatly if you want by running the pipeline yourself (see proteomics_simulator_full_pipline.py).

Basically, it create a data-base of all the possible peptides that are created from a proteomics procedure with proteins defined from the mRNA coding sequences in the 
`<input_mRNAs>` fasta file and their edited versions inferred from the editing_sites_file.

In this command peptides that have more than 3 missed cleavage sites or are shorter than 6 amino-acids or have a mass larger than 5200Da will be discarded.
The pipline dose consider losss of cleavage sites or creation of new cleavage sites due to editing of the coding sequence (inferred from `<editing_sites_file>`)
The dufault enzyme the is used for the in-sillico cleavage is trypsin.

The DB fields names are self-explenatory

After creating the data base, you might want to run another process in order to determine for each of the peptides produced, if it has a uniqe protein source in terms of a gene (and not only isophorm).
This is important because one genes can have many variants (depending on alterative splicing) and ig you want to inffer the expression of an editing site later on, you need to discriminate between peptides that have one gene source (but many different transctipts sources) and peptides that have many gene sources and are not informative when it comes to inferring the presence of an edited (or native) version of a gene in the proteom.

Use this command for that:

`python post_processing_genomic_locations.py <peptides_DB_dir_full_path> <name_of_pepties_file> <path_editing_sites_file>`


You can also create a fasta file containing all the translated proteins (and their edited versions, in all combinations of editing) from the `<input_mRNAs>` from the previous command. This will become handy if you would like to use MaxQuant or other software for the detection of peptides and edited peptides in real proteomics samples (see next stage)

For creating this, you can use the following command:

## Invoking Maxquant for many proteomics samples

In this stage, you will invoke maxquant for many proteomics samples (you should have the .raw files from the Mass-spec machine) simultaniously

Prerequsite for this stage are:

1. MaxQuant software (version 1.6.10.43) installed on your machine

2. Proteome file - fasta file containing the proteom against which you are searching the .raw files for peptides. see previous stage for the creation of such proteom.

Use this command to invoke MaxQuant aginst several proteoms files against your samples .raw files.

`python run_maxquant_sequentialy_for_different_proteomes.py -raw_path <path to samples raw files dir> -proteoms <list of proteoms paths to run MQ against, space seperated> -quantification lfq -o <output suffixess for each of the proteoms resutls>`

use `--help` for more optoins.

This command will invoke MaxQnat for the raw files in `<path to samples raw files dir>` with label free quantification and other default parameters.
MaxQuant results are for each of the proteoms will be stored in `<path to samples raw files dir>` in a folder named "combined_" + suffix_for_proteom 

## Analyzing MQ results

After analyzing several proteomics results, you might want to compare them to the peptides DB created in the simulation stage and quantify the expression of the 


you should have a table of the samples for which you executed MQ for. the tables fields (seperated by tabs) are:
* sample name
* sample sub-name (this is usefull in case you have many samples from one project)
* quantification method (lfq/itraq/tmt)
* phosphoproteomics (True if this is a phosphoproteomics sample, otherwise False)
* sample's tissue
* path to .raw files directory
* notes 1
* notes 2

use this command to analyze MQ results from a certain proteom, stored in the sub-folder `<results_folder>` within the .raw files path:

`python sequential_maxquant_analyses.py -samples_list <path to samples data table> -peptides <path to peptides DB> -sites <path to editing sites table> -combined_output_dir_name <results_folder> -threads 10 -o <suffix for analysis output table>`



