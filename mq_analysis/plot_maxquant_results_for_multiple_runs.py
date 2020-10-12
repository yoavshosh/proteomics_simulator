import re
import os
import sys
import pandas as pd
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statsmodels.stats.multitest as p_adjust
from collections import Counter
from pylab import text
from scipy import stats
from collections import deque
from functools import reduce
from matplotlib import colors
from matplotlib.colors import LogNorm
from heapq import nsmallest
from Bio import SeqIO
import xlsxwriter
from Bio.Alphabet import generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import seaborn as sns
import scipy
from sequential_maxquant_analyses import read_editing_sites_wrt_coding_seqs



all_mm = ['AG','AC','AT','CA','CG','CT','GA','GC','GT','TA','TG','TC']


def get_sites_list(df):
    
    sites = []
    for sample in list(set(df['sample_name'])):
        samples_df = df[df['sample_name'] == sample]
        row_gene_sites = []
        for i, row in samples_df.iterrows():
            for mm_type in eval(row['genomic_keys']):
                row_gene_sites = row_gene_sites + [row['gene_name']+';'+s for s in mm_type]
        sites = sites + list(set(row_gene_sites))
    return sites


def check_possible_contaminant(row):
    contaminants = ['>P00761 SWISS-PROT:P00761|TRYP_PIG Trypsin - Sus scrofa (Pig).',
                     '>Q32MB2 TREMBL:Q32MB2;Q86Y46 Tax_Id=9606 Gene_Symbol=KRT73 Keratin-73',
                     '>P19013 SWISS-PROT:P19013 Tax_Id=9606 Gene_Symbol=KRT4 keratin 4',
                     '>Q7RTT2 TREMBL:Q7RTT2 Tax_Id=9606 Gene_Symbol=KRT78 Keratin-78',
                     '>P15636 SWISS-PROT:P15636 Protease I precursor Lysyl endopeptidase Achromobacter lyticus.',
                     '>P09870 SWISS-PROT:P09870 Arg-C (Clostripain) - Clostridium histolyticum.',
                     '>Q9R4J5 SWISS-PROT:Q9R4J5 Endoproteinase Asp-N',
                     '>P0C1U8 SWISS-PROT:P0C1U8 Endoproteinase Glu-C (V8 -Glutamyl endopeptidase - Staphylococcus aureus.',
                     '>P00766 SWISS-PROT:P00766 Chymotrypsinogen A - Bos taurus (Bovine).',
                     '>P13717 SWISS-PROT:P13717 Nuclease - Serratia marcescens.',
                     '>Q9U6Y5 SWISS-PROT:Q9U6Y5 Green fluorescent protein (GFP-Cter-HisTag)',
                     '>P21578 SWISS-PROT:P21578|LUXY_VIBFI Yellow fluorescent protein (YFP)- Vibrio fischeri.',
                     '>O76009 SWISS-PROT:O76009 Keratin, type I cuticular HA3-I (Hair keratin, type I HA3-I)',
                     '>O76011 SWISS-PROT:O76011 Keratin, type I cuticular HA4 (Hair keratin, type I HA4)',
                     '>O76013 SWISS-PROT:O76013 Keratin, type I cuticular HA6 (Hair keratin, type I HA6)',
                     '>O76014 SWISS-PROT:O76014 Keratin, type I cuticular HA7 (Hair keratin, type I HA7)',
                     '>O76015 SWISS-PROT:O76015 Keratin, type I cuticular HA8 (Hair keratin, type I HA8)',
                     '>P08779 SWISS-PROT:P08779 Tax_Id=9606 Gene_Symbol=KRT16 Keratin, type I cytoskeletal 16',
                     '>Q14525 SWISS-PROT:Q14525 Keratin, type I cuticular HA3-II (Hair keratin, type I HA3-II)',
                     '>Q14532 SWISS-PROT:Q14532 Keratin, type I cuticular HA2 (Hair keratin, type I HA2)',
                     '>Q15323 SWISS-PROT:Q15323 Keratin, type I cuticular HA1 (Hair keratin, type I HA1)',
                     '>Q92764 SWISS-PROT:Q92764 Keratin, type I cuticular HA5 (Hair keratin, type I HA5)',
                     '>Q14533 SWISS-PROT:Q14533 Keratin, type II cuticular Hb1 (Hair keratin, type II Hb1) (ghHKb1) (ghHb1) (MLN 137) - Homo sapiens (Human).',
                     '>Q9NSB4 SWISS-PROT:Q9NSB4 Keratin, type II cuticular Hb2 (Hair keratin, type II Hb2) - Homo sapiens (Human).',
                     '>P78385 SWISS-PROT:P78385 Keratin, type II cuticular Hb3 (Hair keratin, type II Hb3) - Homo sapiens (Human).',
                     '>Q9NSB2 SWISS-PROT:Q9NSB2 Keratin, type II cuticular Hb4 (Hair keratin, type II Hb4) - Homo sapiens (Human).',
                     '>P78386 SWISS-PROT:P78386 Keratin, type II cuticular Hb5 (Hair keratin, type II Hb5) - Homo sapiens (Human).',
                     '>O43790 SWISS-PROT:O43790 Keratin, type II cuticular Hb6 (Hair keratin, type II Hb6) (ghHb6) - Homo sapiens (Human).',
                     '>Q6IFU5 TREMBL:Q6IFU5 Type I hair keratin KA36 - Homo sapiens (Human).',
                     '>Q9UE12 TREMBL:Q9UE12 Type I hair keratin 1 - Homo sapiens (Human).',
                     '>Q8IUT8 TREMBL:Q8IUT8 Type I hair keratin 4 - Homo sapiens (Human).',
                     '>Q6NT21 TREMBL:Q6NT21 Keratin, hair, basic, 3 - Homo sapiens (Human).',
                     '>Q6ISB0 TREMBL:Q6ISB0 Keratin, hair, basic, 4 - Homo sapiens (Human).',
                     '>Q6NTB9 TREMBL:Q6NTB9 Type I hair keratin 3A - Homo sapiens (Human).',
                     '>Q6IFU6 TREMBL:Q6IFU6 Type I hair keratin KA35 - Homo sapiens (Human).',
                     '>P04264 SWISS-PROT:P04264 Tax_Id=9606 Gene_Symbol=KRT1 Keratin, type II cytoskeletal 1',
                     '>P13647 SWISS-PROT:P13647 Tax_Id=9606 Gene_Symbol=KRT5 Keratin, type II cytoskeletal 5',
                     '>P35908 SWISS-PROT:P35908 Tax_Id=9606 Gene_Symbol=KRT2 Keratin, type II cytoskeletal 2 epidermal',
                     '>P13645 SWISS-PROT:P13645 Tax_Id=9606 Gene_Symbol=KRT10 Keratin, type I cytoskeletal 10',
                     '>P35527 SWISS-PROT:P35527 Tax_Id=9606 Gene_Symbol=KRT9 Keratin, type I cytoskeletal 9',
                     '>A3EZ79 TREMBL:A3EZ79;Q32W65;Q6UXC7 Tax_Id=9606 Gene_Symbol=DMKN Dermokine beta-1',
                     '>P02533 SWISS-PROT:P02533 Tax_Id=9606 Gene_Symbol=KRT14 Keratin, type I cytoskeletal 14',
                     '>P02538 SWISS-PROT:P02538 Tax_Id=9606 Gene_Symbol=KRT6A Keratin, type II cytoskeletal 6A',
                     '>P48668 SWISS-PROT:P48668 Tax_Id=9606 Gene_Symbol=KRT6C Keratin, type II cytoskeletal 6C',
                     '>P04259 SWISS-PROT:P04259 Tax_Id=9606 Gene_Symbol=KRT6B Keratin, type II cytoskeletal 6B',
                     '>A3EZ82 TREMBL:A3EZ82;Q32W66 Tax_Id=9606 Gene_Symbol=DMKN Dermokine',
                     '>Q2KIG3 TREMBL:Q2KIG3 (Bos taurus) Similar to carboxypeptidase B2',
                     '>Q0VCM5 TREMBL:Q0VCM5 (Bos taurus) Similar to Inter-alpha-trypsin inhibitor heavy chain H1',
                     '>Q3SZ57 SWISS-PROT:Q3SZ57 (Bos taurus) Alpha-fetoprotein precursor',
                     '>Q9N2I2 SWISS-PROT:Q9N2I2 (Bos taurus) Plasma serine protease inhibitor precursor',
                     '>Q3SZH5 TREMBL:Q3SZH5 (Bos taurus) Similar to Angiotensinogen',
                     '>P28800 SWISS-PROT:P28800 (Bos taurus) Alpha-2-antiplasmin precursor',
                     '>Q1A7A4 TREMBL:Q1A7A4 (Bos taurus) similar to complement component C5',
                     '>P41361 SWISS-PROT:P41361 (Bos taurus) Antithrombin-III precursor',
                     '>Q2YDI2 SWISS-PROT:Q2YDI2 (Bos taurus) Origin recognition complex subunit 4',
                     '>Q3Y5Z3 SWISS-PROT:Q3Y5Z3 (Bos taurus) Adiponectin precursor',
                     '>P81644 SWISS-PROT:P81644 (Bos taurus) Apolipoprotein A-II precursor',
                     '>Q2KJ83 TREMBL:Q2KJ83 (Bos taurus) Similar to Carboxypeptidase N catalytic chain',
                     '>Q2KIT0 TREMBL:Q2KIT0 (Bos taurus) Similar to collagen, type X, alpha 1',
                     '>A2I7N3 TREMBL:A2I7N3;Q27984 (Bos taurus) SERPINA3-7',
                     '>Q3SZV7 TREMBL:Q3SZV7 (Bos taurus) Similar to hemopexin',
                     '>Q2KJC7 TREMBL:Q2KJC7;Q8HZM3 (Bos taurus) Periostin, osteoblast specific factor',
                     '>Q3SZR3 SWISS-PROT:Q3SZR3 (Bos taurus) Alpha-1-acid glycoprotein precursor',
                     '>Q28107 SWISS-PROT:Q28107 (Bos taurus) Coagulation factor V precursor',
                     '>P02672 SWISS-PROT:P02672 (Bos taurus) Fibrinogen alpha chain precursor',
                     '>Q1RMN8 TREMBL:Q1RMN8 (Bos taurus) Similar to Immunoglobulin lambda-like polypeptide 1',
                     '>Q58D62 SWISS-PROT:Q58D62 (Bos taurus) Fetuin-B precursor',
                     '>P06868 SWISS-PROT:P06868 (Bos taurus) Plasminogen precursor',
                     '>Q2KJF1 TREMBL:Q2KJF1 (Bos taurus) Alpha-1-B glycoprotein',
                     '>P02584 SWISS-PROT:P02584 (Bos taurus) Profilin-1',
                     '>P02777 SWISS-PROT:P02777 (Bos taurus) similar to Platelet factor 4',
                     '>Q3SX14 TREMBL:Q3SX14 (Bos taurus) Similar to Gelsolin',
                     '>P17697 SWISS-PROT:P17697 (Bos taurus) Clusterin precursor',
                     '>Q6T181 TREMBL:Q6T181;Q6T182 (Bos taurus) similar to sex hormone-binding globulin',
                     '>P34955 SWISS-PROT:P34955 (Bos taurus) Alpha-1-antiproteinase precursor',
                     '>P21752 SWISS-PROT:P21752 (Bos taurus) Thymosin beta-9',
                     '>Q32PJ2 SWISS-PROT:Q32PJ2 (Bos taurus) Apolipoprotein A-IV precursor',
                     '>Q28194 TREMBL:Q28194 (Bos taurus) Thrombospondin-1',
                     '>P00978 SWISS-PROT:P00978 (Bos taurus) AMBP protein precursor',
                     '>Q5XQN5 SWISS-PROT:Q5XQN5 (Bos taurus) Keratin, type II cytoskeletal 5',
                     '>Q32PI4 TREMBL:Q32PI4 (Bos taurus) Similar to complement factor I',
                     '>Q9TTE1 SWISS-PROT:Q9TTE1 (Bos taurus) Endopin-1 precursor',
                     '>Q2KIU3 TREMBL:Q2KIU3 (Bos taurus) Similar to C1q and tumor necrosis factor related protein 5',
                     '>P01044-1 SWISS-PROT:P01044-1 (Bos taurus) Isoform HMW of Kininogen-1 precursor',
                     '>P67983 SWISS-PROT:P67983 (Bos taurus) Metallothionein-1A',
                     '>Q28065 SWISS-PROT:Q28065 (Bos taurus) C4b-binding protein alpha chain precursor',
                     '>Q862S4 TREMBL:Q862S4 (Bos taurus) Similar to pro alpha 1(I) collagen (Fragment)',
                     '>Q2KIF2 TREMBL:Q2KIF2 (Bos taurus) Similar to leucine-rich alpha-2-glycoprotein 1',
                     '>Q3SX28 TREMBL:Q3SX28;Q5KR48 (Bos taurus) Tropomyosin 2',
                     '>Q0V8M9 TREMBL:Q0V8M9;Q9TRI0 (Bos taurus) similar to inter-alpha (globulin) inhibitor H3 isoform 2',
                     '>Q148H6 TREMBL:Q148H6 (Bos taurus) Hypothetical protein MGC139876',
                     '>Q29RQ1 SWISS-PROT:Q29RQ1 (Bos taurus) Complement component C7 precursor',
                     '>Q95M17 SWISS-PROT:Q95M17 (Bos taurus) Acidic mammalian chitinase precursor',
                     '>P07224 SWISS-PROT:P07224 (Bos taurus) Vitamin K-dependent protein S precursor',
                     '>Q2HJF0 TREMBL:Q2HJF0 (Bos taurus) Similar to Serotransferrin',
                     '>Q2KIH2 TREMBL:Q2KIH2;Q68RU0 (Bos taurus) Ovarian and testicular apolipoprotein N',
                     '>P13646-1 SWISS-PROT:P13646-1 Tax_Id=9606 Gene_Symbol=KRT13 Isoform 1 of Keratin, type I cytoskeletal 13',
                     '>Q04695 SWISS-PROT:Q04695 Tax_Id=9606 Gene_Symbol=KRT17 Keratin, type I cytoskeletal 17',
                     '>A2I7N0 TREMBL:A2I7N0;Q28922;Q3ZEJ6 (Bos taurus) SERPINA3-4',
                     '>P12763 SWISS-PROT:P12763 (Bos taurus) Alpha-2-HS-glycoprotein precursor',
                     '>P17690 SWISS-PROT:P17690 (Bos taurus) Beta-2-glycoprotein 1 precursor',
                     '>P02769 SWISS-PROT:P02769 (Bos taurus) Bovine serum albumin precursor',
                     '>P02676 SWISS-PROT:P02676 (Bos taurus) similar to Fibrinogen beta chain precursor',
                     '>P50448 SWISS-PROT:P50448 (Bos taurus) Factor XIIa inhibitor precursor',
                     '>P01030 SWISS-PROT:P01030 (Bos taurus) similar to Complement C4-A precursor',
                     '>P01966 SWISS-PROT:P01966 (Bos taurus) Hemoglobin subunit alpha',
                     '>P02768-1 SWISS-PROT:P02768-1 Tax_Id=9606 Gene_Symbol=ALB Isoform 1 of Serum albumin precursor',
                     '>P00735 SWISS-PROT:P00735 (Bos taurus) Prothrombin precursor (Fragment)',
                     '>Q03247 SWISS-PROT:Q03247 (Bos taurus) Apolipoprotein E precursor',
                     '>Q3ZBS7 TREMBL:Q3ZBS7 (Bos taurus) Vitronectin',
                     '>Q2UVX4 SWISS-PROT:Q2UVX4 (Bos taurus) Complement C3 precursor',
                     '>Q9TT36 SWISS-PROT:Q9TT36 (Bos taurus) Thyroxine-binding globulin precursor',
                     '>Q28085 SWISS-PROT:Q28085 (Bos taurus) Complement factor H precursor',
                     '>Q3SX09 TREMBL:Q3SX09 (Bos taurus) similar to HBG protein',
                     '>P01045-1 SWISS-PROT:P01045-1 (Bos taurus) Isoform HMW of Kininogen-2 precursor',
                     '>Q3ZBD7 SWISS-PROT:Q3ZBD7 (Bos taurus) Glucose-6-phosphate isomerase',
                     '>Q3MHN2 SWISS-PROT:Q3MHN2 (Bos taurus) Complement component C9 precursor',
                     '>Q9TRI1 TREMBL:Q9TRI1 (Bos taurus) similar to inter-alpha-trypsin inhibitor heavy chain2',
                     '>P15497 SWISS-PROT:P15497 (Bos taurus) Apolipoprotein A-I precursor',
                     '>Q95121 SWISS-PROT:Q95121 (Bos taurus) Pigment epithelium-derived factor precursor',
                     '>Q05443 SWISS-PROT:Q05443 (Bos taurus) Lumican precursor',
                     '>P02070 SWISS-PROT:P02070 (Bos taurus) Hemoglobin subunit beta',
                     '>Q2KIS7 SWISS-PROT:Q2KIS7 (Bos taurus) Tetranectin precursor',
                     '>Q3MHH8 TREMBL:Q3MHH8 (Bos taurus) Amylase, alpha 2B; pancreatic',
                     '>Q3T052 TREMBL:Q3T052;Q5EA67 (Bos taurus) Inter-alpha (Globulin) inhibitor H4',
                     '>Q3KUS7 TREMBL:Q3KUS7 (Bos taurus) Complement factor B',
                     '>Q1RMK2 TREMBL:Q1RMK2 (Bos taurus) IGHM protein',
                     '>Q2TBQ1 TREMBL:Q2TBQ1 (Bos taurus) Coagulation factor XIII, B polypeptide',
                     '>Q05B55 TREMBL:Q05B55 (Bos taurus) Similar to Ig kappa chain C region',
                     '>A2I7N1 TREMBL:A2I7N1 (Bos taurus) SERPINA3-5',
                     '>P04258 SWISS-PROT:P04258 (Bos taurus) Similar to Collagen alpha 1(III) chain',
                     '>Q2KJ62 TREMBL:Q2KJ62 (Bos taurus) KNG protein',
                     '>Q0IIK2 TREMBL:Q0IIK2 (Bos taurus) Transferrin',
                     '>Q3MHN5 SWISS-PROT:Q3MHN5 (Bos taurus) Vitamin D-binding protein precursor',
                     '>P02662 SWISS-PROT:P02662 Alpha-S1-casein - Bos taurus (Bovine).',
                     '>P02663 SWISS-PROT:P02663 Alpha-S2-casein [Contains: Casocidin-1 - Bos taurus (Bovine).',
                     '>P02666 SWISS-PROT:P02666 Beta-casein - Bos taurus (Bovine).',
                     '>P02668 SWISS-PROT:P02668 Kappa-casein [Contains: Casoxin C; Casoxin 6; Casoxin A; Casoxin B; Casoplatelin] - Bos taurus (Bovine).',
                     '>P31096 SWISS-PROT:P31096 Osteopontin - Bos taurus (Bovine).',
                     '>P02754 SWISS-PROT:P02754 Beta-lactoglobulin - Bos taurus (Bovine).',
                     '>P00711 SWISS-PROT:P00711 Alpha-lactalbumin - Bos taurus (Bovine).',
                     '>P62894 SWISS-PROT:P62894 Cytochrome c - Bos taurus (Bovine).',
                     '>Q29443 SWISS-PROT:Q29443 Serotransferrin - Bos taurus (Bovine).',
                     '>P19001 SWISS-PROT:P19001 Tax_Id=10090 Gene_Symbol=Krt19 Keratin, type I cytoskeletal 19',
                     '>A2AB72 TREMBL:A2AB72 Tax_Id=10090 Gene_Symbol=Krt32 Keratin complex 1, acidic, gene 2',
                     '>Q8VED5 TREMBL:Q8VED5 Tax_Id=10090 Gene_Symbol=Krt79 Keratin 79',
                     '>Q61726 TREMBL:Q61726 Tax_Id=10090 Gene_Symbol=5430421N21Rik Keratin type II',
                     '>Q3ZAW8 TREMBL:Q3ZAW8;Q9EQD6;Q9EQD7 Tax_Id=10090 Gene_Symbol=Krt16 Keratin intermediate filament 16a',
                     '>P50446 SWISS-PROT:P50446 Tax_Id=10090 Gene_Symbol=Krt6a Keratin, type II cytoskeletal 6A',
                     '>Q497I4 TREMBL:Q497I4 Tax_Id=10090 Gene_Symbol=Krt35 Krt35 protein',
                     '>Q9D312 SWISS-PROT:Q9D312 Tax_Id=10090 Gene_Symbol=Krt20 Keratin, type I cytoskeletal 20',
                     '>P08730-1 SWISS-PROT:P08730-1 Tax_Id=10090 Gene_Symbol=Krt13 Isoform 1 of Keratin, type I cytoskeletal 13',
                     '>Q922U2 SWISS-PROT:Q922U2 Tax_Id=10090 Gene_Symbol=Krt5 Keratin, type II cytoskeletal 5',
                     '>Q8BGZ7 TREMBL:Q8BGZ7;Q99MH7 Tax_Id=10090 Gene_Symbol=Krt75 10 days neonate skin cDNA, RIKEN full-length enriched library, clone:4732475I03 product:CYTOKERATIN homolog',
                     '>A2A4G1 TREMBL:A2A4G1 Tax_Id=10090 Gene_Symbol=Krt15 Keratin complex 1, acidic, gene 15',
                     '>Q9QWL7 SWISS-PROT:Q9QWL7 Tax_Id=10090 Gene_Symbol=Krt17 Keratin, type I cytoskeletal 17',
                     '>Q6IME9 TREMBL:Q6IME9 Tax_Id=10090 Gene_Symbol=Krt72 Type-II keratin Kb35',
                     '>Q6NXH9 TREMBL:Q6NXH9 Tax_Id=10090 Gene_Symbol=Krt73 Keratin 73',
                     '>A2VCT4 TREMBL:A2VCT4;Q0VDS4;Q6IFX4 Tax_Id=10090 Gene_Symbol=Krt39 keratin 39',
                     '>P07744 SWISS-PROT:P07744 Tax_Id=10090 Gene_Symbol=Krt4 Keratin, type II cytoskeletal 4',
                     '>Q6IFZ6 SWISS-PROT:Q6IFZ6 Tax_Id=10090 Gene_Symbol=Krt77 Keratin, type II cytoskeletal 1b',
                     '>Q6IFX2 SWISS-PROT:Q6IFX2 Tax_Id=10090 Gene_Symbol=Krt42 Keratin, type I cytoskeletal 42',
                     '>Q9R0H5 SWISS-PROT:Q9R0H5 Tax_Id=10090 Gene_Symbol=Krt71 Keratin, type II cytoskeletal 6G',
                     '>Q3TTY5 SWISS-PROT:Q3TTY5 Tax_Id=10090 Gene_Symbol=Krt2 Keratin, type II cytoskeletal 2 epidermal',
                     '>Q0VBK2 TREMBL:Q0VBK2;Q8C1M7 Tax_Id=10090 Gene_Symbol=Krt80 Keratin 80',
                     '>P02535-1 SWISS-PROT:P02535-1 Tax_Id=10090 Gene_Symbol=Krt10 Isoform 1 of Keratin, type I cytoskeletal 10',
                     '>Q61782 TREMBL:Q61782 Tax_Id=10090 Gene_Symbol=Krt14 Type I epidermal keratin (Fragment)',
                     '>A2A5Y0 TREMBL:A2A5Y0 Tax_Id=10090 Gene_Symbol=Krt31 keratin complex 1, acidic, gene 1',
                     '>Q99PS0 SWISS-PROT:Q99PS0 Tax_Id=10090 Gene_Symbol=Krt23 Keratin, type I cytoskeletal 23',
                     '>Q9D646 TREMBL:Q9D646 Tax_Id=10090 Gene_Symbol=Krt34 10 days neonate skin cDNA, RIKEN full-length enriched library, clone:4733401E01 product:keratin complex 1, acidic, gene 4, full insert sequence',
                     '>P05784 SWISS-PROT:P05784 Tax_Id=10090 Gene_Symbol=Krt18 Keratin, type I cytoskeletal 18',
                     '>Q9DCV7 SWISS-PROT:Q9DCV7 Tax_Id=10090 Gene_Symbol=Krt7 Keratin, type II cytoskeletal 7',
                     '>Q9Z2K1 SWISS-PROT:Q9Z2K1 Tax_Id=10090 Gene_Symbol=Krt16 Keratin, type I cytoskeletal 16',
                     '>P07477 SWISS-PROT:P07477 Tax_Id=9606 Gene_Symbol=PRSS1 Trypsin-1 precursor',
                     '>P05787 SWISS-PROT:P05787 Tax_Id=9606 Gene_Symbol=KRT8 Keratin, type II cytoskeletal 8',
                     '>Q6KB66-1 SWISS-PROT:Q6KB66-1 Tax_Id=9606 Gene_Symbol=KRT80 Isoform 1 of Keratin, type II cytoskeletal 80',
                     '>Q7Z794 SWISS-PROT:Q7Z794 Tax_Id=9606 Gene_Symbol=KRT77 Keratin 77',
                     '>Q9BYR9 SWISS-PROT:Q9BYR9 Tax_Id=9606 Gene_Symbol=KRTAP2-4;LOC644350;LOC728934;KAP2.1B;LOC730755;LOC728285 Keratin-associated protein 2-4',
                     '>Q9BYQ5 SWISS-PROT:Q9BYQ5 Tax_Id=9606 Gene_Symbol=KRTAP4-6 Keratin-associated protein 4-6',
                     '>Q9BYR8 SWISS-PROT:Q9BYR8 Tax_Id=9606 Gene_Symbol=KRTAP3-1;LOC100132802 Keratin-associated protein 3-1',
                     '>Q9BYQ7 SWISS-PROT:Q9BYQ7 Tax_Id=9606 Gene_Symbol=KRTAP4-1 Uncharacterized protein ENSP00000381489',
                     '>Q3LI72 SWISS-PROT:Q3LI72 Tax_Id=9606 Gene_Symbol=KRTAP19-5 Keratin-associated protein 19-5',
                     '>Q9BYR4 SWISS-PROT:Q9BYR4 Tax_Id=9606 Gene_Symbol=KRTAP4-3 Keratin-associated protein 4-3',
                     '>Q9BYQ8 SWISS-PROT:Q9BYQ8 Tax_Id=9606 Gene_Symbol=KRTAP4-9 Keratin-associated protein 4-9',
                     '>P60413 SWISS-PROT:P60413 Tax_Id=9606 Gene_Symbol=KRTAP10-12 Keratin-associated protein 10-12',
                     '>P19012 SWISS-PROT:P19012 Tax_Id=9606 Gene_Symbol=KRT15 Keratin, type I cytoskeletal 15',
                     '>Q2M2I5 SWISS-PROT:Q2M2I5 Tax_Id=9606 Gene_Symbol=KRT24 Keratin, type I cytoskeletal 24',
                     '>O95678 SWISS-PROT:O95678 Tax_Id=9606 Gene_Symbol=KRT75 Keratin, type II cytoskeletal 75',
                     '>Q01546 SWISS-PROT:Q01546 Tax_Id=9606 Gene_Symbol=KRT76 Keratin, type II cytoskeletal 2 oral',
                     '>Q99456 SWISS-PROT:Q99456 Tax_Id=9606 Gene_Symbol=KRT12 Keratin, type I cytoskeletal 12',
                     '>Q9H552 TREMBL:Q9H552 Tax_Id=9606 Gene_Symbol=- Keratin-8-like protein 1',
                     '>P35900 SWISS-PROT:P35900 Tax_Id=9606 Gene_Symbol=KRT20 Keratin, type I cytoskeletal 20',
                     '>Q3SY84 SWISS-PROT:Q3SY84 Tax_Id=9606 Gene_Symbol=KRT71 Keratin, type II cytoskeletal 71',
                     '>Q8N1A0 TREMBL:Q8N1A0 Tax_Id=9606 Gene_Symbol=KRT222P truncated type I keratin KA21',
                     '>Q8N1N4-2 SWISS-PROT:Q8N1N4-2 Tax_Id=9606 Gene_Symbol=KRT78 Isoform 2 of Keratin, type II cytoskeletal 78',
                     '>Q5XKE5 SWISS-PROT:Q5XKE5 Tax_Id=9606 Gene_Symbol=KRT79 Keratin, type II cytoskeletal 79',
                     '>P12035 SWISS-PROT:P12035 Tax_Id=9606 Gene_Symbol=KRT3 Keratin, type II cytoskeletal 3',
                     '>Q9C075 SWISS-PROT:Q9C075 Tax_Id=9606 Gene_Symbol=KRT23 Keratin, type I cytoskeletal 23',
                     '>P08729 SWISS-PROT:P08729 Tax_Id=9606 Gene_Symbol=KRT7 Keratin, type II cytoskeletal 7',
                     '>Q7Z3Y8 SWISS-PROT:Q7Z3Y8 Tax_Id=9606 Gene_Symbol=KRT27 Keratin, type I cytoskeletal 27',
                     '>Q7RTS7 SWISS-PROT:Q7RTS7 Tax_Id=9606 Gene_Symbol=KRT74 Keratin, type II cytoskeletal 74',
                     '>Q7Z3Y9 SWISS-PROT:Q7Z3Y9 Tax_Id=9606 Gene_Symbol=KRT26 Keratin, type I cytoskeletal 26',
                     '>Q7Z3Z0 SWISS-PROT:Q7Z3Z0 Tax_Id=9606 Gene_Symbol=KRT25 Keratin, type I cytoskeletal 25',
                     '>Q7Z3Y7 SWISS-PROT:Q7Z3Y7 Tax_Id=9606 Gene_Symbol=KRT28 Keratin 25D',
                     '>P08727 SWISS-PROT:P08727 Tax_Id=9606 Gene_Symbol=KRT19 Keratin, type I cytoskeletal 19',
                     '>Q14CN4-1 SWISS-PROT:Q14CN4-1 Tax_Id=9606 Gene_Symbol=KRT72 Isoform 1 of Keratin, type II cytoskeletal 72',
                     '>Q3KNV1 TREMBL:Q3KNV1;Q96GE1 Tax_Id=9606 Gene_Symbol=KRT7 keratin 7',
                     '>Q86YZ3 SWISS-PROT:Q86YZ3 Tax_Id=9606 Gene_Symbol=HRNR Hornerin',
                     '>P20930 SWISS-PROT:P20930 Tax_Id=9606 Gene_Symbol=FLG Filaggrin',
                     '>Q5D862 SWISS-PROT:Q5D862 Tax_Id=9606 Gene_Symbol=FLG2 Filaggrin-2',
                     '>Streptavidin (S.avidinii)',
                     '>REFSEQ:XP_986630 Tax_Id=10090 Gene_Symbol=Krt33b keratin complex 1, acidic, gene 3',
                     '>REFSEQ:XP_001474382 Tax_Id=10090 Gene_Symbol=1110025L11Rik RIKEN cDNA 1110025L11',
                     '>REFSEQ:XP_092267 Tax_Id=9606 Gene_Symbol=LOC150739 similar to Keratin, type II cytoskeletal 8',
                     '>REFSEQ:XP_932229 Tax_Id=9606 Gene_Symbol=KRT126P similar to Keratin, type II cytoskeletal 2 oral',
                     '>H-INV:HIT000016045 Tax_Id=9606 Gene_Symbol=- Similar to Keratin, type II cytoskeletal 8',
                     '>H-INV:HIT000292931 Tax_Id=9606 Gene_Symbol=- Similar to Keratin, type II cytoskeletal 8',
                     '>H-INV:HIT000015463 Tax_Id=9606 Gene_Symbol=PTPN14 Similar to Keratin 18',
                     '>ENSEMBL:ENSP00000377550 Tax_Id=9606 Gene_Symbol=KRT13 46 kDa protein',
                     '>ENSEMBL:ENSBTAP00000006074 (Bos taurus) 81 kDa protein',
                     '>ENSEMBL:ENSBTAP00000038329 (Bos taurus) 9 kDa protein',
                     '>REFSEQ:XP_001252647 (Bos taurus) similar to endopin 2B',
                     '>ENSEMBL:ENSBTAP00000007350 (Bos taurus) similar to Complement C4-A precursor',
                     '>ENSEMBL:ENSBTAP00000038253 (Bos taurus) 63 kDa protein',
                     '>ENSEMBL:ENSBTAP00000023402 (Bos taurus) 46 kDa protein',
                     '>ENSEMBL:ENSBTAP00000024466 (Bos taurus) 44 kDa protein',
                     '>ENSEMBL:ENSBTAP00000023055 (Bos taurus) 68 kDa protein',
                     '>ENSEMBL:ENSBTAP00000018229 (Bos taurus) 54 kDa protein',
                     '>ENSEMBL:ENSBTAP00000016046 (Bos taurus) similar to fibulin-1 C isoform 1',
                     '>ENSEMBL:ENSBTAP00000024462 (Bos taurus) 47 kDa protein',
                     '>ENSEMBL:ENSBTAP00000014147 (Bos taurus) 12 kDa protein',
                     '>ENSEMBL:ENSBTAP00000033053 (Bos taurus) 15 kDa protein',
                     '>ENSEMBL:ENSBTAP00000001528 (Bos taurus) similar to intersectin long isoform 4',
                     '>ENSEMBL:ENSBTAP00000037665 (Bos taurus) similar to Pregnancy zone protein, partial',
                     '>ENSEMBL:ENSBTAP00000031900 (Bos taurus) 121 kDa protein',
                     '>ENSEMBL:ENSBTAP00000031360 (Bos taurus) 55 kDa protein',
                     '>ENSEMBL:ENSBTAP00000018574 (Bos taurus) 55 kDa protein',
                     '>ENSEMBL:ENSBTAP00000032840 (Bos taurus) similar to apolipoprotein B, partial',
                     '>ENSEMBL:ENSBTAP00000011227 (Bos taurus) 15 kDa protein',
                     '>ENSEMBL:ENSBTAP00000025008 (Bos taurus) hypothetical protein',
                     '>ENSEMBL:ENSBTAP00000034412 (Bos taurus) similar to C4b-binding protein alpha chain',
                     '>ENSEMBL:ENSBTAP00000013050 (Bos taurus) hypothetical protein',
                     '>ENSEMBL:ENSBTAP00000016285 (Bos taurus) similar to peptidoglycan recognition protein L',
                     '>ENSEMBL:ENSBTAP00000024146 (Bos taurus) similar to alpha-2-macroglobulin isoform 1',
                     '>REFSEQ:XP_585019 (Bos taurus) similar to afamin']
    
    possible_contaminants = 0
    for prot in eval(row['proteins_list']):
        for cont in contaminants:
            if prot in cont:
                possible_contaminants+=1
    return possible_contaminants


def calc_mass_change(aa_change):
    mass_before = ProteinAnalysis(aa_change[0]).molecular_weight()
    mass_after = ProteinAnalysis(aa_change[1]).molecular_weight()
    return float((mass_after-mass_before)/mass_before)


def plot_editings_sites_repetitions(fig_path,names, sites, groups = [1,2,3,4,5,6,7,8,9,10]):
    col = ['b','g','salmon','y','m','tan','r','silver','yellow','crimson','teal','k','orange','brown','limegreen','navy','gold','c','olive']
    
    bars = []
    sites_counts = [list(Counter(s).values()) for s in sites]
    for s in sites_counts:
        b=[]
        for g in sorted(groups[:-1]):
            b.append(len([c for c in s if c==g]))
        b.append(len([c for c in s if c>=groups[-1]]))
        bars.append(b)
    
    # set width of bar
    barWidth = 0.2
    
    groups_labels = []
    for i in groups:
        if i < groups[-1]:
            groups_labels.append(str(i))
    groups_labels.append(str(groups[-1])+'<=')
    

    # Set position of bar on X axis
    rs = [np.arange(1,len(bars[0])+1)]
    r_temp = rs[0]
    for i in names[1:]:
        r_temp = [x + barWidth for x in r_temp]
        rs.append(np.array(r_temp))

     
    # Make the plot
    for i in range(len(bars)):
        plt.bar(list(rs[i]), bars[i], color=col[i], width=barWidth, edgecolor='white', label=names[i])
      
    # Add xticks on the middle of the group bars
    plt.xlabel('Repetitions', fontweight='bold')
    plt.ylabel('Sites', fontweight='bold')
    plt.xticks([g+0.2 for g in groups], groups_labels) 
    plt.yticks(np.arange(0,25,5))
    plt.title('Discovered Editing Sites')
    
    # Create legend & Show graphic
    plt.legend()
    plt.savefig(fig_path)
#    plt.show()
    plt.close()


def scatter_discovered_peptides(df, name, path):
    color_markers = [('y', '+'), ('r', '*'), ('b', '>'), ('g', '<'), ('m', '^'), ('tan', 'D'), ('silver', '.'), ('yellow', 'o'), ('crimson', '<'), ('teal', '>'), ('k', 'D'), ('orange', '^'), ('brown', 'D'), ('limegreen', 'D'), ('navy', '^'), ('gold', 'o'), ('c', '.'), ('olive', '.')]
    
    tissues = Counter(list(df['tissue']))
    
    for i,t in enumerate(tissues.items()):
        print( t[0]+' '+str(t[1]))
        tissues_dfs = df[df['tissue']==t[0]]
        plt.scatter(list(tissues_dfs['total_peptides']),list(tissues_dfs['edited_peptides']), label = t[0]+' '+str(t[1]), color = color_markers[i][0], marker = color_markers[i][1])
    plt.legend(loc='center left', bbox_to_anchor=(0.8, 0.5))
    plt.title('Peptides for Samples - ' + name + ' proteom')
    plt.xlabel('total peptides')
    plt.savefig(path)
    plt.close()
    
   

def plot_editing_level_wrt_detections(final_combined_detections_df,path):


    detection_df = final_combined_detections_df[final_combined_detections_df['proteom_editing_type']=='AG']
    scipy.stats.pearsonr
#    g = sns.jointplot(x = 'edited_detections', y = 'editing_level', data = detection_df, marker='o', color='red')
#    sns_plot = sns.lmplot(x= 'edited_detections', y="editing_level", hue="name", data=detection_df, palette="Set1");
#    sns_plot = sns.jointplot("edited_detections", "editing_level", data=detection_df, kind = 'reg', marginal_kws=dict(bins=15),annot_kws=dict(stat="r"), s=40, edgecolor="w", linewidth=1)
#    sns_plot.savefig(path+'edited_detections_lin_reg.png')

    
    g = sns.lmplot('edited_detections', 'editing_level', data=detection_df, sharex=False, sharey=False, fit_reg = True)
#    ax = g.axes[0,0]
#    g = sns.JointGrid(x="edited_detections", y="editing_level", data=detection_df)
#    g = sns.JointGrid(x="edited_detections", y="editing_level", data=detection_df)
#    g = g.plot_marginals(sns.distplot, kde=False, color=".5")  
#    g = g.plot(plt.scatter, sns.distplot)
#    g = g.plot(sns.lmplot('edited_detections', 'editing_level', data=detection_df, sharex=False, sharey=False, fit_reg = True))
#    g = g.plot(sns.regplot('edited_detections', 'editing_level', data=detection_df))
    g.set(ylim=(0, None))
    g = g.savefig(path+'edited_detections_lin_reg.jpg')
    
    
def plot_edited_native_and_both_detection_wrt_editing_level(path, detection_df, sites_df, minimal_detections, editing_levels=[0.2, 0.4, 0.6, 0.8]):
    
    bars_edited_only = []
    bars_native_only = []
    bars_both = []
    
    if 1 not in editing_levels:
        editing_levels.append(1)
    
    lower_b = 0
    for i,el in enumerate(editing_levels):
        el_detections = detection_df[np.logical_and(detection_df['editing_level']>lower_b,detection_df['editing_level']<=el)]
        both = len(el_detections[el_detections['both_detections']>minimal_detections])
        edited_only = len(el_detections[np.logical_and(el_detections['edited_detections']>=minimal_detections,el_detections['native_detections']==0)])
        native_only = len(el_detections[np.logical_and(el_detections['edited_detections']==0,el_detections['native_detections']>minimal_detections)])
        total_sites = len(list(set(list(sites_df[np.logical_and(sites_df['editing_level']>lower_b,sites_df['editing_level']<el)]['genomic_key_base1']))))
        bars_edited_only.append(edited_only/total_sites)
        bars_native_only.append(native_only/total_sites)
        bars_both.append(both/total_sites)
        lower_b = el
        
    barWidth = 0.25
    r1 = np.arange(len(bars_edited_only))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    
    plt.bar(r1, bars_edited_only, color='b', width=barWidth, edgecolor='white', label='edited')
    plt.bar(r2, bars_native_only, color='r', width=barWidth, edgecolor='white', label='native')
    plt.bar(r3, bars_both, color='g', width=barWidth, edgecolor='white', label='both')
    
    plt.xlabel('editing levels', fontweight='bold')
    plt.ylabel('percent detection', fontweight='bold')
    ticks = ['0-'+str(editing_levels[0])]
    for i,el in enumerate(editing_levels[:-1]):
        ticks.append(str(el)+'-'+str(editing_levels[i+1]))
    plt.xticks([r + barWidth for r in range(len(bars_edited_only))], ticks)
    
    
    plt.legend(bbox_to_anchor=(0.9, 1))
    plt.savefig(path+'detections_per_editing_levels.jpg')
    plt.close()



def plot_editing_level_detections_grouped_bars_all_proteomes(path,final_combined_detections_df,proteoms = ['AG','non_AG','random_AG'],editing_levels=[[0,0.02],[0.02,0.05],[0.05,1]],score_cutoff=None,PEP_cutoff=None):
    
    pval_results={}
    detection_dict = {}
    for p in proteoms:
        print(p)
        df = final_combined_detections_df[final_combined_detections_df['proteom_editing_type']==p]
        if p=='random_AG':
            sites_edited=df[df['total_edited_detections']>0]
            sites_covered=df
            if score_cutoff is not None:
                sites_edited=sites_edited[sites_edited['max_edited_score']>=score_cutoff]
                sites_covered=sites_covered[np.logical_or(sites_covered['max_edited_score']>=score_cutoff,sites_covered['max_non_edited_score']>=score_cutoff)]
            if PEP_cutoff is not None:
                sites_edited=sites_edited[sites_edited['min_edited_PEP']<=PEP_cutoff]
                sites_covered=sites_covered[np.logical_or(sites_covered['min_edited_PEP']<=PEP_cutoff,sites_covered['min_non_edited_PEP']<=PEP_cutoff)]

            detection_dict.update({p:(len(sites_edited),len(sites_covered))})
            print(len(sites_edited))
            print(len(sites_covered))
        else:
            for el in editing_levels:
                print(str(el))
                df_for_el = df[np.logical_and(df['editing_level']>el[0], df['editing_level']<=el[1])]
                sites_edited=df_for_el[df_for_el['total_edited_detections']>0]
                sites_covered=df_for_el
                if score_cutoff is not None:
                    sites_edited=sites_edited[sites_edited['max_edited_score']>score_cutoff]
                    sites_covered=sites_covered[np.logical_or(sites_covered['max_edited_score']>=score_cutoff,sites_covered['max_non_edited_score']>=score_cutoff)]
                if PEP_cutoff is not None:
                    sites_edited=sites_edited[sites_edited['min_edited_PEP']<PEP_cutoff]
                    sites_covered=sites_covered[np.logical_or(sites_covered['min_edited_PEP']<=PEP_cutoff,sites_covered['min_non_edited_PEP']<=PEP_cutoff)]
                    
                if p in detection_dict.keys():
                    detection_dict[p].update({str(el[0])+'-'+str(el[1]):(len(sites_edited),len(sites_covered))})
                else:
                    detection_dict.update({p:{str(el[0])+'-'+str(el[1]):(len(sites_edited),len(sites_covered))}})
                print(len(sites_edited))
                print(len(sites_covered))

    p_val_dict = {}
    for el in editing_levels:
        print(str(el))
        for p in proteoms:
            print(p)
            if p!='random_AG':
                p_val = stats.fisher_exact([[detection_dict['random_AG'][1],detection_dict['random_AG'][0]],[detection_dict[p][str(el[0])+'-'+str(el[1])][1],detection_dict[p][str(el[0])+'-'+str(el[1])][0]]])[1]
                print(str(p_val))
                if p in p_val_dict.keys():
                    p_val_dict[p].update({str(el[0])+'-'+str(el[1]):p_val})
                else:
                    p_val_dict.update({p:{str(el[0])+'-'+str(el[1]):p_val}})
                pval_results.update({p+' '+str(el[0])+'-'+str(el[1]):p_val})
                


    barWidth = 0.25
    bars1 = [detection_dict['AG']['0-0.02'][0]/detection_dict['AG']['0-0.02'][1],detection_dict['non_AG']['0-0.02'][0]/detection_dict['non_AG']['0-0.02'][1], 0]
    bars2 = [detection_dict['AG']['0.02-0.05'][0]/detection_dict['AG']['0.02-0.05'][1],detection_dict['non_AG']['0.02-0.05'][0]/detection_dict['non_AG']['0.02-0.05'][1], detection_dict['random_AG'][0]/detection_dict['random_AG'][1]]
    bars3 = [detection_dict['AG']['0.05-1'][0]/detection_dict['AG']['0.05-1'][1],detection_dict['non_AG']['0.05-1'][0]/detection_dict['non_AG']['0.05-1'][1], 0]
    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    # Make the plot
    b1 = plt.bar(r1, bars1, color='lightskyblue', width=barWidth, edgecolor='white', label='0-0.02')
    b2 = plt.bar(r2, bars2, color=['steelblue', 'steelblue', 'silver'], width=barWidth, edgecolor='white', label='0.02-0.05')
    b3 = plt.bar(r3, bars3, color='navy', width=barWidth, edgecolor='white', label='0.05-1')
     
    # Add xticks on the middle of the group bars
    plt.xlabel('Proteom', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['AG', 'non_AG', 'random_AG'])
    plt.ylabel('Edited Detection Rate', fontweight='bold') 
    
    p_val_texts=[]
    for el in editing_levels:
        for p in proteoms:
            if p=='random_AG':
                p_val_texts.append('')
            else:
                pval=p_val_dict[p][str(el[0])+'-'+str(el[1])]
                if pval>=0.05:
                    p_val_texts.append('ns')
                elif pval<0.0001:
                    p_val_texts.append('****')
                elif pval<0.001:
                    p_val_texts.append('***')
                elif pval<0.01:
                    p_val_texts.append('**')
                elif pval<0.05:
                    p_val_texts.append('*')

    # Create legend & Show graphic
    plt.legend(title='Editing Levels')
    
    max_height=0
    for i,rect in enumerate(b1 + b2 + b3):
        height = rect.get_height()
        if height>max_height:
            max_height=height
        plt.text(rect.get_x() + rect.get_width()/2.0, height, p_val_texts[i], ha='center', va='bottom')
    
    plt.ylim(0,max_height+0.1)
    plt.savefig(path+'Edited_Detection_Rate_PEP_lower_than_'+str(PEP_cutoff)+'.png')
    plt.close()
    
    return pval_results, detection_dict
 
def plot_scores(scores_distributions_dict,path,name,log_scale=False):

    colors=['b','r','g']
    j=0
    for k,v in scores_distributions_dict.items():
        x = np.arange(0,len(v))
        y=list(v)
        plt.plot(x,y,color=colors[j],label=k)
        j+=1
    if log_scale:
        plt.yscale('log')
    plt.legend()
    
    plt.savefig(path+name+'.png')
    plt.close()

def plot_sites_scores_distribution_by_proteoms(path,scores_distributions_dict,name,log_scale=False,hist=False):
    
    sns.distplot(scores_distributions_dict['AG'], color='blue', hist=hist, norm_hist=True, label='AG')
    sns.distplot(scores_distributions_dict['non_AG'], color='red', hist=hist, norm_hist=True, label='non_AG')
    sns.distplot(scores_distributions_dict['random_AG'], color='green', hist=hist, norm_hist=True ,label='random_AG')
    plt.savefig(path+name+'_dist.png')
    plt.close()
  
    
    
def fined_maximal_score_for_sites_covered(row, sites_dfs_dict):
    proteom=row['proteom_editing_type']
    sites_df=sites_dfs_dict[proteom]
    sites_df=sites_df[sites_df['genomic_site']==row['site']]
    edited_sites_dfs=sites_df[sites_df['site_is_edited']]
    non_edited_sites_dfs=sites_df[~sites_df['site_is_edited']]
    if row['total_edited_detections']>0:
        max_edited_score=max(edited_sites_dfs['Score'])
    else:
        max_edited_score=0
    if row['native_detections']>0:
        max_non_edited_score=max(non_edited_sites_dfs['Score'])
    else:
        max_non_edited_score=0
    row['max_edited_score']=max_edited_score
    row['max_non_edited_score']=max_non_edited_score
    return row

def fined_minimal_PEP_for_sites_covered(row, sites_dfs_dict):
    
    proteom=row['proteom_editing_type']
    sites_df=sites_dfs_dict[proteom]
    sites_df=sites_df[sites_df['genomic_site']==row['site']]
    edited_sites_dfs=sites_df[sites_df['site_is_edited']]
    non_edited_sites_dfs=sites_df[~sites_df['site_is_edited']]
    if row['total_edited_detections']>0:
        min_edited_PEP=min(edited_sites_dfs['PEP'])
    else:
        min_edited_PEP=1
    if row['native_detections']>0:
        min_non_edited_PEP=min(non_edited_sites_dfs['PEP'])
    else:
        min_non_edited_PEP=1
    row['min_edited_PEP']=min_edited_PEP
    row['min_non_edited_PEP']=min_non_edited_PEP
    return row

def plot_pvals_results_per_cutoff(path,pvals_per_PEP_dict, cutoffs):
    
    colors = ['deepskyblue','blue','navy','coral','red','darkred']
    i=0
    for k, v in pvals_per_PEP_dict.items():
        plt.plot(cutoffs, v, label=k,color=colors[i])
        i+=1
    plt.legend()
    plt.savefig(path+'pvals_per_PEP.png')
    plt.close()
  

def plot_detections_results_per_cutoff(path,detections_per_PEP_dict, cutoffs):
    
    colors = ['deepskyblue','blue','navy','coral','red','darkred']
    i=0
    for k, v in detections_per_PEP_dict.items():
        edited_detections = [i[0] for i in v]
        plt.plot(cutoffs, edited_detections, label=k,color=colors[i])
        i+=1
    plt.legend()
    plt.savefig(path+'edited_detections_per_PEP.png')
    plt.close()
    
    i=0
    for k, v in detections_per_PEP_dict.items():
        edited_coverage_rate = [i[0]/float(i[1]) for i in v]
        plt.plot(cutoffs, edited_coverage_rate, label=k,color=colors[i])
        i+=1
    plt.legend()
    plt.savefig(path+'edited_coverage_rate_per_PEP.png')
    plt.close()
    
    
if __name__ == '__main__':
        
    # path = 'E:/RNA_editing_Large_files/MQ_human/only_pride/'
    path = 'C:/Users/shosh/Google_Drive/RNA_Editing/Human_editom_expression/results/'
    samples_list_name='organized_samples'
    # sites_file = 'C:/Users/shosh/Google_Drive/RNA_Editing/Human_editom_expression/human_recoding_editom_wrt_to_coding_sequences_with_genomic_coor_with_additional_data.txt'
    
    
    names = ['AG','non_AG','random_AG']
    summ_files = [path+'results_from_'+samples_list_name+'_'+n.lower()+'_editings' for n in names]
    sites_files = [path+'sites_detection_from_'+samples_list_name+'_'+n.lower()+'_editings' for n in names]
    
    quantification_methods = ['lfq','lfq+','unknown','tmt','itraq']
    minimal_detections = 1
    filter_by_quantification = False
    filter_contaminants = True
    filter_sites_by_mass_change = None    #mass increase/decrease due to editing
    score_cutoff=None
    PEP_cutoff=None

    # all_sites_df = read_editing_sites_wrt_coding_seqs(sites_file)
    # all_sites_df['editing_level'] = pd.to_numeric(all_sites_df['editing_level'])
    
    #name for combined analysis
    analysis_name  = ''
    for n in names:
        analysis_name = analysis_name + '_' + n
    if filter_by_quantification:
        analysis_name = analysis_name[1:] + '_LFQ'
    else:
        analysis_name = analysis_name[1:] + '_all_quantifications'
        
#    sys.stdout = open(path + analysis_name + '.txt', 'w')

    #read results summaries per sites-type and store in a dictionary of dataframes
    summaries_dfs_dict = {}
    writer = pd.ExcelWriter(path + 'sammaries.xlsx', engine='xlsxwriter')
    for i,file in enumerate(summ_files):
        df = pd.read_csv(file,sep = '\t', header = 0)
        df['full_sample_name'] = df.apply(lambda row: row['sample_name']+'_'+str(row['sub_name']), axis = 1)
        df.to_excel(writer, sheet_name = names[i], index = False)
        if filter_by_quantification:
            df = df[df['quantification'].isin(quantification_methods)]
        summaries_dfs_dict.update({names[i]:df})
    writer.save()
   
    #read sites dataframes per sites-type and store in a dictionary of dataframes
    sites_dfs_dict = {}
    writer = pd.ExcelWriter(path + 'peptides.xlsx', engine='xlsxwriter')
    for i,file in enumerate(sites_files):
        df = pd.read_csv(file,sep = '\t', header = 0)
        for j, mm in enumerate(all_mm):
            df[mm+'sites'] = df.apply(lambda row: eval(row['genomic_keys'])[j], axis = 1)
        df['full_sample_name'] = df.apply(lambda row: row['sample_name']+'_'+str(row['sub_name']), axis = 1)
        df['mass_change'] = df.apply(lambda row: calc_mass_change(row['aa_change']), axis = 1)
        df.to_excel(writer, sheet_name = names[i], index = False)
        if filter_by_quantification:
            df = df[df['quantification'].isin(quantification_methods)]
        if filter_sites_by_mass_change is not None:
            df = df[np.fabs(df['mass_change']) >= filter_sites_by_mass_change]
        if filter_contaminants:
            df['possible_contaminant']=df.apply(lambda row: check_possible_contaminant(row),axis=1)
            df = df[df['possible_contaminant']==0]
#        if score_cutoff is not None:
#            df = df[df['Score']>=score_cutoff]
        sites_dfs_dict.update({names[i]:df})
    writer.save()
    
    #in all analysis types, consider only results from samples for which all types of analysis were preformed
    finshed_run_in_all_analyses = list(set.intersection(*map(set,[[str(i)+'_'+str(j) for i,j in zip(list(summaries_dfs_dict[n]['sample_name']),list(summaries_dfs_dict[n]['sub_name']))] for n in names])))
    for k,v in summaries_dfs_dict.items():
        summaries_dfs_dict[k] = v[v['full_sample_name'].isin(finshed_run_in_all_analyses)]
    for k,v in sites_dfs_dict.items():
        sites_dfs_dict[k] = v[v['full_sample_name'].isin(finshed_run_in_all_analyses)]
    
    #for each analysis type, print a list of all sites and the number of samples in which edited/native/both versions was found (note that edited and native columns contain counts of sampls for which both versions were found)
    detections_dfs = []
    for k, v in sites_dfs_dict.items():
        name = k
        df = v
        edited_sites = df['genomic_site']
        detections = []
        for site in set(edited_sites):
            site_df = df[df['genomic_site']==site]
            samples_edited = list(set(list(site_df[site_df['site_is_edited']]['full_sample_name'])))
            samples_native = list(set(list(site_df[~site_df['site_is_edited']]['full_sample_name'])))
            edited_and_native = list(set(samples_native).intersection(samples_edited))
            editing_level = list(site_df['editing_level'])[0]
            mm_type = site_df.iloc[0]['mm_type']
            tissues = site_df[site_df['site_is_edited']].groupby('tissue')['full_sample_name'].apply(list)
            
            #calculate detections per tissue
            tissues_detections = []
            for t in set(list(tissues.index)):
                tissues_detections.append(t+'_'+str(len(set(tissues[t]))))    
        
            detections.append((site,mm_type,len(samples_edited),len(samples_native),len(edited_and_native),len(samples_edited)+len(edited_and_native),editing_level,tissues_detections))
        
        detections.sort(key = lambda x: (x[5], x[2], x[3], x[4],len(x[7])), reverse = True)
        detection_df = pd.DataFrame(data = detections, columns = ['site','mm_type','edited_detections','native_detections','both_detections','total_edited_detections','editing_level','tissues_detection'])
        
        edited_detections = len(detection_df[detection_df['edited_detections']>0])
        print('\n\nSites Edited and Native Versions Detections for ' + name + ' Analysis (' + str(edited_detections)+ ' edited detections)' )
        print(detection_df.to_string(index = False, justify = 'left'))
        
        detection_df['proteom_editing_type'] = detection_df.apply(lambda x: k, axis = 1)
        detections_dfs.append(detection_df)
        detection_df[np.logical_or(detection_df['edited_detections']>=minimal_detections,detection_df['native_detections']>=minimal_detections)]
    
    maximal_scores_distributions_dict={}
    for n in ['AG','non_AG','random_AG']:
        df=sites_dfs_dict[n]
        df=df[df['site_is_edited']]
        df['max_score']=df.apply(lambda row: max(df[df['genomic_site']==row['genomic_site']]['Score']),axis=1)
        df=df.sort_values('max_score',ascending='False')
        df=df.drop_duplicates('genomic_site')
        maximal_scores_distributions_dict.update({n:df['max_score']})
    plot_sites_scores_distribution_by_proteoms(path,maximal_scores_distributions_dict,'score',log_scale=False)  
    plot_scores(maximal_scores_distributions_dict,path,'score',log_scale=False)
    
    minimal_PEP_distributions_dict={}
    for n in ['AG','non_AG','random_AG']:
        df=sites_dfs_dict[n]
        df=df[df['site_is_edited']]
        df['min_PEP']=df.apply(lambda row: min(df[df['genomic_site']==row['genomic_site']]['PEP']),axis=1)
        df=df.sort_values('min_PEP',ascending='True')
        df=df.drop_duplicates('genomic_site')
        minimal_PEP_distributions_dict.update({n:df['min_PEP']})
    plot_sites_scores_distribution_by_proteoms(path,minimal_PEP_distributions_dict,'PEP',log_scale=True)  
    plot_scores(minimal_PEP_distributions_dict,path,'PEP',log_scale=True)
    
    
    final_combined_detections_df = pd.concat(detections_dfs)    
    final_combined_detections_df=final_combined_detections_df.apply(lambda row: fined_maximal_score_for_sites_covered(row,sites_dfs_dict), axis=1)
    final_combined_detections_df=final_combined_detections_df.apply(lambda row: fined_minimal_PEP_for_sites_covered(row,sites_dfs_dict), axis=1)
    writer = pd.ExcelWriter(path + 'supported_sites_all_proteoms.xlsx', engine='xlsxwriter')
    final_combined_detections_df.to_excel(writer, sheet_name = 'sites', index = False)
    writer.save()
  
    PEP_edited_random_ag = final_combined_detections_df[final_combined_detections_df['proteom_editing_type']=='random_AG']
    PEP_edited_random_ag = PEP_edited_random_ag[PEP_edited_random_ag['edited_detections']>0]
    
    PEP_results_columns = ['PEP_upper_bound','proteom','editing_level','edited_detection','all_detection','pval']
    data = []
    # PEP_values = set(sorted(PEP_edited_random_ag['min_edited_PEP'].values))
    PEP_values = [0.02]
    
    pvals_per_PEP_dict = {}
    detections_per_PEP_dict = {}
    editing_levels_for_calc = [[0,0.02],[0.02,0.05],[0.05,1]]
    proteoms = ['AG','non_AG','random_AG']
    for p in proteoms:    
        for el in editing_levels_for_calc:
            if p!='random_AG':
                detections_per_PEP_dict.update({p+' '+str(el[0])+'-'+str(el[1]):[]})
                pvals_per_PEP_dict.update({p+' '+str(el[0])+'-'+str(el[1]):[]})
    for cut in PEP_values:
        pval_results, detection_dict = detections_results = plot_editing_level_detections_grouped_bars_all_proteomes(path,final_combined_detections_df,proteoms = proteoms,editing_levels=editing_levels_for_calc,score_cutoff=None,PEP_cutoff=cut)
        for p in proteoms:    
            if p!='random_AG':
                for el in editing_levels_for_calc:
                    detections_per_PEP_dict[p+' '+str(el[0])+'-'+str(el[1])] += [detection_dict[p][str(el[0])+'-'+str(el[1])]]
                    pvals_per_PEP_dict[p+' '+str(el[0])+'-'+str(el[1])] += [pval_results[p+' '+str(el[0])+'-'+str(el[1])]]
                    data.append((cut,p,str(el[0])+'-'+str(el[1]),detection_dict[p][str(el[0])+'-'+str(el[1])][0],detection_dict[p][str(el[0])+'-'+str(el[1])][1],pval_results[p+' '+str(el[0])+'-'+str(el[1])],))
            else:
                data.append((cut,p,'-',detection_dict[p][0],detection_dict[p][1],'-',))

    results_per_PEP_df = pd.DataFrame(data = data, columns=PEP_results_columns)
    results_per_PEP_df.to_excel(path+'results_per_PEP.xlsx',index=False)
    
    
    score_results_columns = ['score_lower_bound','proteom','editing_level','edited_detection','all_detection','pval']
    data = []
    score_values = set(sorted(PEP_edited_random_ag['max_edited_score'].values))
    pvals_per_score_dict = {}
    detections_per_score_dict = {}
    editing_levels_for_calc = [[0,0.02],[0.02,0.05],[0.05,1]]
    proteoms = ['AG','non_AG','random_AG']
    for p in proteoms:    
        for el in editing_levels_for_calc:
            if p!='random_AG':
                detections_per_score_dict.update({p+' '+str(el[0])+'-'+str(el[1]):[]})
                pvals_per_score_dict.update({p+' '+str(el[0])+'-'+str(el[1]):[]})
    for cut in score_values:
        pval_results, detection_dict = detections_results = plot_editing_level_detections_grouped_bars_all_proteomes(path,final_combined_detections_df,proteoms = proteoms,editing_levels=editing_levels_for_calc,PEP_cutoff=None,score_cutoff=cut)
        for p in proteoms:    
            if p!='random_AG':
                for el in editing_levels_for_calc:    
                    detections_per_score_dict[p+' '+str(el[0])+'-'+str(el[1])] += [detection_dict[p][str(el[0])+'-'+str(el[1])]]
                    pvals_per_score_dict[p+' '+str(el[0])+'-'+str(el[1])] += [pval_results[p+' '+str(el[0])+'-'+str(el[1])]]
                    data.append((cut,p,str(el[0])+'-'+str(el[1]),detection_dict[p][str(el[0])+'-'+str(el[1])][0],detection_dict[p][str(el[0])+'-'+str(el[1])][1],pval_results[p+' '+str(el[0])+'-'+str(el[1])],))
            else:
                data.append((cut,p,'-',detection_dict[p][0],detection_dict[p][1],'-',))

    results_per_score_df = pd.DataFrame(data = data, columns=score_results_columns)
    results_per_score_df.to_excel(path+'results_per_score.xlsx',index=False)
    

#    plot_edited_native_and_both_detection_wrt_editing_level(path,final_combined_detections_df[final_combined_detections_df['proteom_editing_type']=='AG'],all_sites_df,minimal_detections,score_cutoff=score_cutoff)
    
    
    #metadata spread sheet of samples finished in all proteoms
    combined_table_dfs=[]
    for n in ['AG','non_AG','random_AG']:
        if n=='AG':
            df=summaries_dfs_dict[n][['sample_name','sub_name','quantification','phosphoproteom','tissue','total_peptides','edited_peptides','AG_detected','total_editing_sitse']]
            df=df.rename({'sample_name':'ID','sub_name':'sample','AG_detected':'AG_edited_sites','total_editing_sitse':'total_edited_AG_sites'}, axis='columns')
            combined_table_dfs.append(df)
        elif n=='non_AG':
            df=summaries_dfs_dict[n][['sample_name','sub_name','total_peptides','edited_peptides','AC_detected','AT_detected','CA_detected','CG_detected','CT_detected','GA_detected','GC_detected','GT_detected','TA_detected','TG_detected','TC_detected','total_editing_sitse']]
            df=df.rename({'sample_name':'ID','sub_name':'sample','total_editing_sitse':'total_edited_non_AG_sites'}, axis='columns')
            for c in df.columns:
                if '_detected' in c:
                    df=df.rename({c:c.replace('detected','edited_sites')}, axis='columns')
            combined_table_dfs.append(df)
        elif n=='random_AG':
            df=summaries_dfs_dict[n][['sample_name','sub_name','total_peptides','edited_peptides','AG_detected','total_editing_sitse']]
            df=df.rename({'sample_name':'ID','sub_name':'sample','AG_detected':'random_AG_edited_sites','total_editing_sitse':'total_edited_random_AG_sites'}, axis='columns')
            combined_table_dfs.append(df)
    
    df1=combined_table_dfs[0].merge(combined_table_dfs[1],on=['ID','sample'],suffixes=('_AG', '_non_AG'))
    df2=df1.merge(combined_table_dfs[2],on=['ID','sample'],suffixes=('', '_random_AG'))
    df2['ID']=df2.apply(lambda row: row['ID'].split('_')[0] if 'PXD' in row['ID'] else row['ID'], axis=1)
    df2['source']=df2.apply(lambda row: 'PRIDE' if 'PXD' in row['ID'] else 'CPTAC', axis=1)
    df2=df2.rename({'total_peptides':'total_peptides_random_AG','edited_peptides':'edited_peptides_random_AG'},axis='columns')
    df2=df2[['source','ID','sample','quantification','phosphoproteom','tissue','total_peptides_AG','edited_peptides_AG','total_edited_AG_sites','total_peptides_non_AG','edited_peptides_non_AG','total_edited_non_AG_sites','total_peptides_random_AG','edited_peptides_random_AG','total_edited_random_AG_sites','AG_edited_sites','random_AG_edited_sites','AC_edited_sites','AT_edited_sites','CA_edited_sites','CG_edited_sites','CT_edited_sites','GA_edited_sites','GC_edited_sites','GT_edited_sites','TA_edited_sites','TG_edited_sites','TC_edited_sites']]
    df2.to_excel(path+'samples_metadata.xlsx') 
    
    
