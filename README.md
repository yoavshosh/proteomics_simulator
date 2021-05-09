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

s
You can also have consensus transcripts from some published genome and coordinates of mutation within that genome for which you would like to asses the protein expression.