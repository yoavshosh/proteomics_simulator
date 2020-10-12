import re
import itertools as it
import sys
from rna_peptide import rna_peptide
from collections import deque
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna, IUPAC, generic_protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#sys.stdout = open('C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/' +'test3.txt', 'w')


#input_path = 'C:/Users/user/Google_Drive/RNA_Editing/proteomics_simulator/test_files/'
##input_fasta = 'sim.fasta'
#input_fasta = 'in_frame_rna_from_orfs_squ.fasta'
digestion_rule = re.compile(r'(CGC|CGA|CGG|CGT|AGA|AGG|AAA|AAG)(?!(CCC|CCA|CCT|CCG))|(?<=TGG)(AAA|AAG)(?=(CCC|CCA|CCT|CCG))|(?<=ATG)(CGC|CGA|CGG|CGT|AGA|AGG)(?=(CCC|CCA|CCT|CCG))')
#sequence = 'AAACGCTAGAAACCACGGAAACCCACTATCAGGCATTCAGGCATCAGCATGCGCCCTATTGCACGC'
#sequence = 'ACTACCACCACTACTACTACTACACAAATTTTCTATGGCCATGTAATCATCATGTAGTAAGTAGTGAGAAGAGATGGTGCGGCAAACGAGGTATCTCCTGGTTGGGAATATACCCGAAAAAATCACCAAGGAGAAAATCACAGAACACTTTAAGAACTATGGAAAGATCCAGAATGTAAAGTTTCATCCAAAAAAAGAAAATGAACTTGGCATTACAGCAACTGTGGCATTTATGGACATTAAATGTGCCTCTAAGGCATATAGTTCAGAGAACAAAATAGAAAATGTAGTGCTTCGGACAGAATACAGTGAGGGCACTGCAACAGGAAGTGTTGTTACAAGACCTGTGGACCCTGGATTGCCTACAACTGGAGCCCACAGGGGACAGTACCTACATGCTCGAGGGCCATCTACTAGCTTCACAACACGGACTAAAGGTAGTCTTGATGAAGATGATACATACAGGGAATATTACTATGGTGGTAGCGGCAGCACTGGCAGTAGGGAAAACAGTCAATTTGAAAACAGATATCCCACATATGATGAACAAGTACAGCATCAGTCAAGGTCTTATCAAAGAAATAGGACCTCATTTTCCAAGCAGTATGAGCGAGGACCACAACGACAGAAGTCTAGCAGTGCGAGTGTCACGGGTGTTGTCACTAGCGGTGGTAGTGGCAGTAGTGTCACAGGGGTTGGGAATAGCAATACAGTTGGCAGTACAGTGACACCCAGCAGTGGTGGTAGTGTGGTTGGGAGTGGAGCTGTTGGAACCAGTGGCAGTGGGGGTGGGGGTGCTGGTAGTGGTGCTGGCAGTGGAGGAGGCAACAATAGTGGAGGTGGTGGCACCAACAATGGAGGGGACACTGGCAACACTGGGGTCAGTGTTGTGGTCTCTAATACCTTGCAAAAGCATTCCTTTGATGAGCGTTACCACGAACAGTATAGCAGCGACAGAGAAAGTTTTGAAAGGCCACGGATATCCACATTGTCTCCAACATGTTCCATCTCTCAGGCCAGCAGCCGACACAGTAGTGGTCGGCGGAAAGGTTCTGTCTCTGTGTCCCGTTCGTCCAGTCGCAGTCGGTCAAGGTCCAGTTCCAGCAGTGAGTATAGTCGGTCTTCATCAGGAAGCAGTGGATCAAGATCACATTCCAGAAGTCGTGCAAGCAGCCGGGGGAGTTTGCGGCGAATACCTGTTGGTTCAGCATCGAGTCAAAAACACTCTAAATCAGAAAGTCAAAGTATGATGTCCAGTAAATCCTCAACAGCAGGATCGCAATCGGAAAGAAAACCACTTGGGATTTGTGTAACGGGTTTACCATTAAGATCCAGTGATACTAGTCTTCGAGATGGTCTTTTTCATGAATACAAAAAACATGGAAAAGTTACTTCTGTACAGGTGCTAGGACAAGGCGAAGAAAGATATTCAGTGGTCAGTTTCCGGAAACCTGAAGATGCGGCAAAGGCATTAGAAGCTTCACAGGACAAAATGTTTTTTGGAAAGAAAATCAAAGTTATGGCTCATGAAGGAGTTGAAGTGGAAGACAATGAGTTCCGTCCCCCAGAAGCTGAATGTGATGAATACCACCCTAAGGCAACACGAACGTTATTTGTGGGAAACTTAGAAAAAGAAATAACGACACAAGAGCTTAGAGATCGTTTTAAGTCTTTTGGAGAAATTATAGATATTGATATCAAACGTCAGGGTGCAGTGTCAGCATATGCATTTGTACAGTATGCAGATATTACAAGTGTTGTAAAATCACTCCGAAAAATGGAGGGCGAACATATAGGAGCCAATAAAATAAAGCTTGGTTTTGGTAAAAGTATGCCAACCAATTGTGTGTGGCTCGATAACCTGTCTGAATCCGTAACAGAGAAGTTTCTTTGCCAACAATTTGGCCGTTATGGCGAGGTCAGTCATGCAGTGATTGATCGGTATAAAGGAAGGGCCTTGATATATTTCACAAGTATGGATACAGCTCAATATGCAGTCACTGAGATGAGAAATAGAATTTTGAAGAAGAAACGAGTGCAGACTGGTGACTTTAGACTCGGTGGACATCTGAATATTGATTTTGCTAGTCGAGATTGTCAAACAGCTTTCTTTGAAAAAATGGAGCGGACGGGTCAACTGCAACCTGGAGAACGTCCTGATGAACGAGGAGGGCGCCTCTACCGGCACAATAGGGGAGCCAATTTCGAGCAATATTATGATTCTCCTGCAACCCCTTCCGCTTACAAGGAAGATAGCGCAAGATACGAGACTGGCACTGTAACACCTGTTGCTACAACTGCACGAGTTGGACTGACCCCACCACAGTTCAGTACGACAAAACGTAATCGAGCAACAAATTTTCAAGCCTCTTCGACCAGACAAGGATCCTCATACAGAGGAAGACATTTTACAGAAAGTACTTACTCAGAAGAATTCCCGGCAAACCAGAGACAGAGACATATTGATGAATACAGTCAAGGAAGTGCTCCTCCTTATGCTGAGGAAGATTCCTATGAACAGGAGCTGAGAGAATATGGCTACAGTCAACGGGAACGCCGTGAGAGGGGTAACGTAACTCCACCGAGACGATACCACAGTCCTTATCGGGAAGATAGGCATTCAAATGCCGACAAAGAATCGACGCGGGGTAGCCGTTATGAAGAAGCCTTCACCAAAGACAGATATGATTATTACAGTAGTGACAGTCCTGTCGTCAGTGACAATGGTCATGATGAGAACTGGGTGGATCATCATACATTGGATCTTTCAGATCGTCATGCCACACCAACAAAAATGAAATATTCAAGTCGGAACGAGACACGCTTAAAGCAATCGTCTCGGGATCACAGCCCTGCCATTGAGTTTTCCTACAGGCACTCACCAGTGAGACGAAGAAGCAAGAGCAGGAGCAGAAGTAAAAGTGGGAGTCGGACAGAAAGTCCTAGTCGAAGTAGGAGCCGGACGCCTCTAAGGGAAAGGCGGAAAGAATCTGTGTTCCGGCCACGTTCTCCCCACACACCACAAACACCTCCTACACCAACAGAAGAAGAACCTCCCATTTATTCTAAAATCCGATCAATCAGCCCCACACCCCGGTTTGAGGAATTTGAGCGAACAGGAGAATTGAAAAAAATATGCAGGACAGATGATAAAATTGAAGGCATTAGGTCTTCGGATCATGAAAAGGTGTACAGTATAAAATGCTTGGGTACTCCATCTATGTCTTTCATGAAAGGGGAAACGACACAAGCCCGTAAAAAGAGTAGGGAGTTTGATGATAGCAAACATGATGTAATAGCAGGACGTCTGCATAAGTCCCAGGTGCCTGCAACGAACCGGCTGTTAGAATCTGCGCTCCAAGACAGACATCAAGTCTTGTTAAAGAAGGCAGATGAAAACCGTCTCAATAGACGCAAACTATTAGACCAGCCACATGGTAGCTCGAGTCATTCCAGTGGCCACGATGATTCAGATCCTGGTTCTTTTGCTGAGACTGATCTTGGCCATTTGCATCGAGAAAAGCGACTACTCCTAGAAAAATTAAAACAACTGGAAGACGCTGGAAGCCCCAGTGATAATGATTCTTTACCGACACATGATGATCTGCGAGAAGATGGACACGGGCGATTTGTACCGAGGAAAGGTCGACCTGAGGATCCTCATTGTAGCAGTTTACGGCATAGTCTGAGCGACATACACGCAAAGGAAAAACGAAAACTAGAAACTTCTCAGAGTTTTCGAAAACAGATGGAAGCCAGGCGCCTGTTAGAACAGACTCAGGGAGAGTCAAATTTGTTATCAAAACCAAAAAAGAGTGTTGAAAAAGTTGTGATGCTTGATGGTATGGATAGCGAAGGAGAGGAGGAAATAAGTATTACCTCTCCGAAAGCAATCCCCACTGCATTCAAGCCAAAGAGAAAGAAGAAACTCGAAGAAGACCCAGTGTCTGGTGTTCGTACTGGCAGGTTTTATCGCATGAAGAAAAGTGGAATGATGAGCAGTGATGACGAAAAAGATCATTCAAAATCTCTCGAAGCAGAAGAACCTGTTATGCGACGCACCAGTCGCGAGAGCAAAGACGAACTATCCGCTGTTATAAGCCGGTTTAAGGAGCAGTCAAGACATTCTGAATTAAGCCGGGAGGTTTCCGAGTCTGGGTCTCACAAGGAATCCAGAGAGTATAAAATGGTCCATATGACCTCGAAGGATGGAATGGAGACGGGTACAAAAATTATCAAAAAAGAACTATCAGTCGACGATGTCGTTCAGGATCCTCGACATGAAATAAAGAAAGAGGAACCTCTTTCTTTACCGTTACCCCGTTTTGCTCTTAAGCACAACACGTCACCTATTGAGTCTCCTCTTCCTTGTCCTTCACCACCTCTATTGCCTGTTGGGCCCAAGGCTCGTTCTCCTTGTTTCTCTCCACAATGCAGTAGCAAATCCCAGTCTCCGGCAATGTCACCTGCTGGTTCGTTATCGGATGTCTCTCATGATACAGTTAATTCTAATCTCCCCGAGTCTGGACTGAACAGGTTACACGAGCAAAAAGATGAAGATCAGCTTAGTAAAGATTTGGCTGCCAGTTCAAGCGGAATCAAAACCCACGTACCAAACCCTGATGGTGACGACAATTTATCAGACAGTTTGAGTGAACTTAGCTGTTCCTCTCTGGATGAAAAAATTCGACAGTTGGATAAAAAACTGAGTATGACTCCTGTTCCACGGCCAGTTGATACAATCACCCCTGGGCTTTATACGAAGTTTAAAATCAAAAAGAAGGAAACACCTCTGTCTGGTGGAAGCATGCTGACATCTAGTCTGCGGAGTGAGCCATCTGATATTGTCAAATCATTGCTTTCTCGGTCGAGTATATTTGATCAAGATTCAAAGCGGTTAGAGCAGATAAATGAAAAATACAAACCAAAAGAGATCAATATAAATATTGAAGGCTCCCCACCCAAAATGAATATCCGTACTAAAGCAGCGGCCAAAGAAATGCCCCCTACAATGCCTCCAAATTTCAGCCCTCTCATGCCCTTTGGTCCTTTTAGTTCCACTTGTAACATTACTTCCCCTTTTCAGCCTTTTACCTCTAACTTGTTCCCGCAGCACCCTGTGATACGCTTGTCAAATACCCCTTCTACTCCACAGTTTTGGGGAGGTTCAACTCCAACCATGAGCTCGCCAACAAAACCACTCGTGGATGTGGCAGCGCCAGTTTCAGTGTTAAAGAAGACTCCGCTCACTCCTGCCACGACAGAGCTGTCGGTCAAAACTGAAGAGTTACCCCATGCTTTGCTTTCTTCGCCACCACCATCTGTGGCAGAAATTAGTCTGCCCCCTGTACGGCAGCCAGTGGATACTCCACTTGCAAGTTCACACACTCCTATGGTGAGTCCATGTTCGGAACAAACTGAACCTGCAGTGAAGAAGGAACCTAGCGTCCCTTTAGTAACTGGTAAAAATGAAGTTGCCGTTTGGATGTCCCCAAGGAAAGATCTCAAAGTAGAAGAGATGAGTGAAGTTGATCCCAAAGAAAAAAGAGAATTTATCCCTGGTCCAGTATCTGTGCCGTCCACATCATCTGGTTCCTGTCTGGGAAAGAGGAAATCCTCAGATGATTTGGATTCCCCCAGGAATAAGATGATGAAAACTGAACCCAAAACTTGTCCTCTTGCAATGAAGCCTCAGGAAGTTCCAAAGAACGACAGCTTTGAGAAGACCGAAGCAAAAGACGAGCCCGTGCAGCAGCCTCGAAAGAAAGCCAAGCACGGCGAGTCTACAAAACCCAAAGAGGAAAAAGAGTGTGAGAAGAAATCTCCTGGGGTACCTTCAGTCAAGCCCAGCAAATCGCCTGAGAAGAAACTGGTAAAGACACCAAGTCGTGAAAAAAACGAAAAAAAGGAGCCAACAGAGGGCGAGAAGAAAGATGTGGAAAAGAAACCAGCTGAACTACCAACTCCGAAGAAAGATCAAAAAGCTTCAAAAGACAGTGATAAAAGCGCTGAGAAATGGAAATCAGGAGGCAAGGAACCAGTCTGTGATTTAGATGAGAAGGCTAAGTTGAAGAATAGCCACTCTGAACGGAAGCATGATAAAAGTGAAAAGAAAGAAAAGAACCGTACAAAAAGCTCAGAGGGGAAGCCTGACAAGGAGAGAAGCTCCTCAAAAGAAGGAGAGGCAAAGGAGTCTCACAAAAAGGACAAATTAAAATCTAAAACTGACGGTCGGTCAAAAACCAAATCGGAAGATAAAGGAGAAGGGGATGCTTCCAAAAATGGGCGGGGTGGCCATGAAAAAATAGAAGCTGAGAATCTATTAAAAGCCTCCAAGTCTCACCGGAAATCTGAAGAAGACAAAGGAGGAGATGACAAAACGGACAATGTAAAATGCTCTGAAGTTAACAAGATTGTTAAATCCCATCATCGTATTGAGCCACCAGAGGGCGAAGCGTCAACAAAAAACAATAAGACACATCACAAAACAGACCCAGAAAAAAGTGAACTTGAAACCTGCAAGGCGGGTAAGTCTCACCATAAGGCTGACAATGATAAAATAGATGGGGAAAGTGTGAAGTCTTCTGGAAAAACTCACCAGAAACCAGCAGAGGGCCTGGAAAAAGTTGAGGCTGATGCCGCCAAGAATAAATGTAACAAAGTGGAAACTGAAAAGACAGAACCAGAAAATACTAAGGCCAGTAAATCTCACCATAAATATGAATGTGAAAAGGGTGAGGCCGAAGGAACAAAACCTAGCAGCAAGTCTCACCACAAACACGACCAATTAAAGTCTGAGGCAGACAGCTTGAAATCGAGTAAGTTGCATCACAAGTCTGACAAAAGTAACCGCAGTTCTTCAGATCGGACGCGTAAAGATAGCGTTAAAAGCAGCTCTGATCGTAAGAAGAACGAAAAAGAGAAAATGGACAAATCTGAAAAACAGGAGAAGTCTGATCGGTTAGAAAAAATGGATAAACCTGACAAAGATACCAAGACCGAAAACAGGGAAAAGTTAGAAAAGACTGAAAAGTCCGAGAAAGAAAAATTTGATAAAGAAAAATCCGATAAAACGGAAAAGGAAAAATCTGACAAGGAAAAGTCCGAGAAAGTGGAAAAGGCGGAAAAAGAAAAGTCCGAAAAGGAAAAATCAGACAAGGAGAAAGATAAAGCTGAAAAGGAGAAGCAAGATAAGAGCAAACCGGATAAAGAAAAACCGGAAAAAGAAAAACCAGATAAGGACAAACAAGAAAAACCCGACAAAGACAAGCAGGAAAAACCTGAGAAAGAGAAACCTGACAAATCTGAAAAGGAAAAACCAGACAAAGAAAAACCCGATCGATCAGAAAAGGATCGGTCAGAAAAGGACAAGTCCGACAGGTCTGAAAAATCTGACAGATCTGAGAAAACTGAGAGCACCGAAAAGATGGACAAGTCCGAGAAAGAAAAGTCCGATAAGGAAAAGTGTGGAGATGAAAAGAAGAGCGAGTTCTGCCTCGATGCCTTTGGTCCGTATGTGTCCATGTACGACAAGGTGAAGAGACGATCCTGCAGTAACAAGGACAAAGACATGGAAGACATGCGCAAGAAACTGAGCCAACTTAAAAACACTCGTAAGAAAAGAGGCAGCAAACCAAGTCGAAGTGATGAGACAGAGAGCAGCGTACAGAACACTGATGATGAGAGTAGTACAAGCAATTCGTTTGTGGCTTCCAAGAAGGAGACGCTAACAGAAGAAAAAGATGAAACTACAAACAAAAAGAAAAAACGAAATGTAATCGAGAGCTCCTCTTCAGAAGAAGACCCTCCAACATTTATACCCACACCCTCCAAGAGAGAACTGAAGAAATCTTCTCCAAAGACCCGACCTATGCACAAAAAGTCCAAGGCCGTCTTAGAAGTCTCTTCTGATAGTGATACTGAAGAGAAAGATTATAAAAAGATCTCATCTAAACTGGCCAAGCAGAAGTTGGGCTCATCATCAACGTTGCTGCATGACCAGAGAAGCAACTCCACCGATGCTAAAAAAGACAGTCAGGATGAAATGATGACAGAAAAGACAGACAAACCGAATCGGAAAAAGTCTAGTAAAGAAAAGAAGAAGAAATGCGAAAAAGAGAAGAAGAAATCCAAGAAGTCTGCTGACAGGTCTTCTAAAGATGAGAAGAAGGTTTTACCAAAATCTTCACCTGTGTTCAATGATATCAGTGAAATGCCAGGCATGGACTCACAAGGTTCCCAGGTTGTTAACATCCTGTCGGACAACATATACATGTCCATAGCCCATGAAGATTCTCCCAAGGAGGATTATCCCCAGATGCAAAAATCAAATAAGCCATCAGTGGCAACAGCATCAACAGTGGAAGACGAAGAAGATGAGGAAGAAAATGCTTCTGATGTTGCTATTGGCACCACTCCTGGCACGAATACCGCAACTGCAGATGAAGACGAGTCCAAACCTTTTGCCCAGCGGAAATCAAGCAAAGATTCCTTCCCTCACAAAATAGAAGAAAGCCATCAGGTTGAGAGCTTGTCAGAGAAACTCTTCCTACCTATAAGTCTGAGTGAGGCTACAGCCAAGGACTCCGAGGAGAACCACCACGCAGATACAGAGGATGAGCCCCCTGACATGTCATCCAAATCGAACAAGAGCAAGAAGAAGAAAAAGAAAGAGAAGAAGAAACATGATACTAAGTCTAAGAAAAGTTCCAAGAAAAACATTGTGAAGATCTCAATGCTCGATGACAATGTTTTCCTAAATGATCAGCCTGAGGTTACGCAAGAGGAAGAGGATGAAGAAGAGGAAGAAGAAGGTACAAAGTTCTTCCAGATTTCTAGCCTTGAGTCTTTTTTTCTTTCTGGCTCAAAAGAAGAACCAAATGAGCTCATAAAAACCGAGCCCCTTGTCGCAGAGGAAGATATGAAAGAAGAAACCGATGTGGGTAGTCTAAAAAGGACGTCATCTTGCGAAGATATCAAGGTAGCAACTAACGATGACTATTTGGAAGAGCCTGTAAAAGAAGAGGCTTCCACAGATTACATGGAGCCTGAAGACCCGCTTCCTACCAAGGTTGAAGAGCCACCAAAGGAAGAACTGGAGGCTGAAGAAGATGAAGTTGAAGAGGAAATCGAGGAGATCAATAAACCAAATACCGACTCCTCCCCGCTGGAAGAAGACAACTCCAACAATTTAGAAAAATGTGAAAGTTTTCAAGATGACGTCTACGATTTTGACAAGTCGGAAAATATCACTGAGAAACTGACTGATTGTGAACAGTTGAAGAAATTTGAAAAAATCTACGACTCTGAGAAGGTGAAAAAGGTTTCTGATGAAGCTGTCGATTCTGACATTCCAAAAAGTGTTGGAGAAGAGTCCCCTGAGTTCCTCTTATCACCAGGAGAGGAGGTCACTGAGGAAGACCCTGATCAACCAGAACCCACACCTGAAGAGCAGTCAGATAACGATGTTCAAAAAGAAACAACAGACGATATTCAAGCTCAAAGCCCTGAGCACTCAGATACAAAAGAGTTGGACAGTCCTTATCAGTTAAAAGAACAAGAAGATACTGTATCGGATGTGTTTACTGATTTGAAGAAACTCGATAAAAGCTGCTCTGAAGTTGAAGACAACAGTGAAGATAAGATTGTTGCAAAAAAGCCTGAAAATGAGGATTTACTTGAAAAAGAAGAACCAGAAAAGACTGAAGAAGAATCGAAATCTCCTGAAACCTCAGAAGTTGAAAAAGATCCGACAGACGAATCTTGTGAGACAACTGTAGATGGTGAGGATCTGAAACGAGATGAAGGCAAAGATCAAAATTCCAGCAGTTTTGAGAAACTAGCCGAAACCAACGAAAAGGTAATTGGTATGTTCAACTTTGGCCAATTACAAGTTCATGTTGTTGACAAGAATACAGGTAACCTTCACAAGGTGCATGACGCCGAGTTATCCGATAAAGATAGTCGGGCCGCTCTTGGTGTGGATCAGGTTCTTGACAAAGCTAACAGCAACACTGATGCCAACAGGAATGTAGACCACATCTATGACTTTGACGATACCACTGAAGAAATGTACGGTAATACTTATCTCGATAAAGGAGATAGCAAATCCAGACAGTTCAACTACAAAAAGGTCTTTATGATGCAACCCCGTAAATCGGACAGAATTTATAATTCTGATGTCACCAAGAAAATCTGCGAGAAGATTGAGCAGAAGATTTACAGCCTAGAAACGCCAACCAAATTTGATACCATCTTCAAGAAAATTAACAAAAAGGCGTACAATATGGACGAGGCTAATTCTACTGAAGATAAGATTTATGAGGATGAGATTTGCAAGACGCGCAATCGAATTGCAGATCTCCACCACCCGCAGCAGCATACAGATAAAACAGATGGTGACGTCTATGATTTCGATGAACAAGACAAGCTTGAAGTTAATGAATTGTTAGATTCCAAGAAAACGGATCATATAATCCATGATCTTCACGATTACAAGAAAACCGAAGAGCATTTACCCAATTTTGATGGAGGCAAGAAAAGTACCGAAAAACTGGCAAGGTTCCAAGATTGGTTTACAAATGACGATAAGATATATTACTTTGACCCTGTTAAAGATTTTGGTGCCAAGCCTAACAAAAGCCGTTCCAAATGTGGTAAACACAAATCTGGTCAGGGGAAACGCCGCAATAAAGGTAAAATTGTTACCTTTACTTGGCCTCGTGAACTAGACTTAAAATCATCGGTTGAATCTGCAACGGACAATGATACTGTCACTGCCGATACGTTTGAGATGAAAGAGAGCGATTATGCTGAACTCGATGACCAAGAAGATAAACTTGTCATCGATACAAGTGACATGTGCGTCTATGATGAAAATTCAACCGAACCGATTGAAAAAGAGGCTGATAAAGCAAAGCTGTGTAGCTCGGAAACAAAAACCGAAAACGATTTAACAGCCAATGATGAAAAAAACGAAGATATCCAAGACGACGTGAAAGAAAAGGTGTCAAGTGAAGAGGAAGAAAAAGAGCAAGAAGAAACGGAGGCTCTGGATACCACGTCTTTTAAAAATAATCAAGACACTTTTGGAATGTCGATCAAACAGAAAGTGGTGGTGAATTATGGACTAACAAATCACAGTGGCACAACCCCCGACAAGATAGAATCGGTTATAGGCTTTGGAATGCTATCTGGCAACAAAGACCTGCTGTCGGCATCTAAATTGTCTGATGACTGTGACAGCAGCGATGAAGGTGCGCTCAGAATCGATGAAGATGAAGTAGACGAAAATACTGATAGGGATGAAATGACACTTGAAGACAAAACCCCAGAAAAAGAGAAAGAAATTGAAGAAGAAGAAAAAGTGATCAAGAAGCCTGAACCCTATGAAGAAGAGGAGGAGGAAGAAGC'
#sequence = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
#sequence = str(rec[0].seq)
#missed_cleavages = 1
#min_length = 6
#max_length = 100
#a2g_sites = [32, 123, 155, 176, 192, 207, 208, 231, 247, 248, 255, 260, 270, 279, 281, 286, 308, 310, 376, 380, 472, 474, 475, 485, 488, 517, 521, 531, 532, 553, 564, 581, 623, 641, 698, 705, 708, 733, 764, 773, 874, 875, 1007, 1008, 1039, 1063, 1070, 1092, 1119, 1206, 1221, 1232, 1235, 1258, 1274, 1280, 1336, 1337, 1353, 1358, 1366, 1369, 1383, 1386, 1446, 1475, 1478, 1480, 1484, 1488, 1489, 1490, 1497, 1498, 1503, 1526, 1544, 1613, 1662, 1717, 1830, 1859, 1877, 1898, 1922, 1925, 1939, 1949, 1952, 2003, 2017, 2033, 2039, 2040, 2054, 2056, 2066, 2069, 2075, 2084, 2099, 2122, 2132, 2141, 2153, 2157, 2160, 2168, 2178, 2188, 2196, 2200, 2208, 2212, 2220, 2226, 2228, 2229, 2237, 2240, 2245, 2253, 2261, 2280, 2313, 2316, 2343, 2346, 2370, 2373, 2387, 2399, 2400, 2414, 2420, 2423, 2483, 2489, 2498, 2586, 2624, 2634, 2642, 2645, 2646, 2717, 2790, 2802, 2822, 2823, 2837, 2886, 2923, 2924, 3036, 3064, 3070, 3087, 3217, 3292, 3335, 3572, 3697, 3700, 3715, 3805, 3810, 3813, 3839, 3892, 3928, 4106, 4126, 4156, 4185, 4470, 4475, 4669, 4670, 5073, 5141, 5142, 5201, 5202, 5226, 5247, 5374, 5513, 5534, 5544, 5584, 5626, 5631, 5688, 5731, 5761, 5788, 5919, 5969, 5974, 5975, 6050, 6051, 6071, 6106, 6108, 6135, 6142, 6147, 6183, 6254, 6284, 6332, 6349, 6352, 6361, 6383, 6384, 6423, 6461, 6464, 6469, 6473, 6503, 6966, 6998, 7008, 7040, 7043, 7061, 7082, 7088, 7229, 7328, 7336, 7340, 7348, 7363, 7421, 7422, 7439, 7460, 7479, 7480, 7481, 7482, 7507, 7546, 7590, 7628, 7643, 7650, 7667, 7757, 7836, 7872, 7877, 8064, 8081, 8091, 8092, 8093, 8094, 8125, 8126, 8128, 8151, 8170, 8176, 8192, 8299, 8502, 8509, 8510, 8517, 8522, 8585, 8636, 8754, 8810, 8987, 9056, 9092, 9096, 9119, 9128, 9131, 9174, 9182]
#a2g_sites = []
#c2t_sites = []
#a2g_sites = [1,30,38]
#c2t_sites = [4,15,55]
#a2g_sites = [2, 7, 11, 16, 18, 20, 24]
#c2t_sites = [3, 12, 13, 22, 25]
#a_to_g_editings = a_list
##c_to_t_editings = c_list
#look_ahead_for_editnig = 3
#look_behind_for_editing = 6

#sys.stdout = open(input_path +'test3.txt', 'w')

def create_edited_rna_peptides(sequence,a2g_sites,c2t_sites,digestion_rule, missed_cleavages = 0,
                               look_ahead_for_editnig=3,look_behind_for_editing=6,min_length = None, max_mass = None, max_sites_per_pep = None):
    
    
    a2g_sites = sorted(a2g_sites)
    c2t_sites = sorted(c2t_sites)
    rna_peptides_list = []
    under_min_length = []
    over_max_sites_peps = []
    over_mass_peps =[]
    over_mass_pep_combs_dict = {}
    
    seq_length = len(sequence)
    cleavage_sites, binary_for_fixed_cs = create_all_cleavage_sites(sequence,digestion_rule,a2g_sites,c2t_sites,look_ahead_for_editnig,look_behind_for_editing)
    cleavage_sites_lengeth = len(cleavage_sites)

    
    #running over all posiblle windows from a site until mc+1 fixed sites ahead
    for i in range(cleavage_sites_lengeth-1):
        fixed_sites_cnt = 0
        j=1
        
        #running over all j (1+j -> index of cleavage site at the end of the subsequence window) untill betwin i-th cs and (i+j)-th cleavage site there are mc number of fixed site or i+j site is sequence end
        #in this method all possible windows for subsequences are examined
        while (fixed_sites_cnt < missed_cleavages+1) and (cleavage_sites_lengeth-i-1>=j):    
            enter_seq = True
            examin_each_comb_mass = False
            edge = None
            append_zero = False
            append_seq_end = False
            subsequence_coor = [cleavage_sites[i],cleavage_sites[i+j]-1]
            coor_str = '_'+str(subsequence_coor[0])+'_'+str(subsequence_coor[1])
            subsequence = sequence[subsequence_coor[0]:subsequence_coor[1]+1]
            
            if binary_for_fixed_cs[i+j]: #count site if site ahead is fixed.
                fixed_sites_cnt+=1
            
            if min_length != None: #rule out peptide if shorter than min_length
                if len(subsequence) < min_length:
                    enter_seq = False
                    under_min_length.append(subsequence_coor)
        
# =============================================================================
            #determining range for editing perutations:
            if binary_for_fixed_cs[i]: #site is fixed and therfore upstream editing dosent matter
                edit_range_start = cleavage_sites[i]
            else:    
                if cleavage_sites[i]-look_behind_for_editing < look_behind_for_editing: #sequence is close to beginning, looking for permutations from beginning
                    edit_range_start = 0
                else:
                    edit_range_start = cleavage_sites[i]-look_behind_for_editing
             
            if binary_for_fixed_cs[i+j]: #site is fixed and therfore downstream editing dosent matter
                edit_range_end = cleavage_sites[i+j]
            else:    
                if cleavage_sites[i+j]+look_ahead_for_editnig > seq_length: #sequence is close to end, looking for permutations till end
                    edit_range_end = seq_length
                else:
                    edit_range_end = cleavage_sites[i+j]+look_ahead_for_editnig
             
            permutation_range = (edit_range_start,edit_range_end-1) #this is the range in wich ro examine editing combinations (could extand further than subsequence (peptide) coordinates)
# =============================================================================
            
            # getting all editing sites in permutations window
            a2g_sites_in_perm_range = find_sites_in_window(permutation_range,a2g_sites)
            c2t_sites_in_perm_range = find_sites_in_window(permutation_range,c2t_sites)
            number_of_sites_in_perm_range = len(a2g_sites_in_perm_range) + len(c2t_sites_in_perm_range)

            #if too many editing sites in peptide, dont examine peptide
            if max_sites_per_pep != None:
                if number_of_sites_in_perm_range > max_sites_per_pep:
                    enter_seq = False
                    over_max_sites_peps.append(subsequence_coor)
                    
                    
            #last condition to rule out peptide - lower bound mass > max_mass allowed
            #getting all editing sites in susequence window (only for getting maximal\minimal mass range)
            minimal_mass, maximal_mass = get_mass_range(subsequence,subsequence_coor,
                                                        find_sites_in_window(subsequence_coor,a2g_sites),
                                                        find_sites_in_window(subsequence_coor,c2t_sites))
            
            if max_mass != None:
                if minimal_mass > max_mass:
                    enter_seq = False
                    over_mass_peps.append(subsequence_coor)
                    #if max_mass is betwin peptide minimal mass and maxiamal mass each editing comb mass is calc before appended
                elif maximal_mass > max_mass: 
                    examin_each_comb_mass = True
# =============================================================================
            if enter_seq:    
                sebsequence_in_perm_range = sequence[edit_range_start:edit_range_end] #on this subsequence we examine different combinations (could extand further than subsequence (peptide) coordinates)
                 
                #all eiditing sites permutation in permutation_range
                permutations = find_permutations_in_range(creat_editing_sites_permutations_one_type(a2g_sites_in_perm_range),
                                                          creat_editing_sites_permutations_one_type(c2t_sites_in_perm_range))
             
                #determining if sequence beginning or end is in cleavage_sites
                if edit_range_start == 0:
                    append_zero = True
                    edge = 5
                if edit_range_end == seq_length:
                    append_seq_end = True
                    edge = 3
                
                for perm in permutations:
                    extand_sebsequence_to_edit = edit_rna_as_peptide(sebsequence_in_perm_range,permutation_range,perm[0],perm[1])
 
                    #the edited cleavage sites withing the range of subsequence (including_bounderies)
                    perm_cs = [k+edit_range_start for k in find_cleavage_sites(digestion_rule,extand_sebsequence_to_edit,append_zero=append_zero, append_seq_end=append_seq_end) 
                                if (k+edit_range_start>=cleavage_sites[i] and k+edit_range_start<=cleavage_sites[i+j])]
                    
                    #if cleavage site at start is fixed, append to cs list for perm
                    if edit_range_start == cleavage_sites[i] and edit_range_start not in perm_cs:
                        perm_cs.insert(0,edit_range_start)
                    #if cleavage site at end is fixed, append to cs list for perm
                    if edit_range_end == cleavage_sites[i+j] and edit_range_end not in perm_cs:
                        perm_cs.append(edit_range_end)
                     
                     
                    #if bouderies still exists as cleavage sites and no more than mc cleavage sites in perm_cs - a valid seq
                    if set([cleavage_sites[i], cleavage_sites[i+j]]).issubset(perm_cs) and len(perm_cs)<=missed_cleavages+2:
                        edited_sub_seq = edit_rna_as_peptide(subsequence,subsequence_coor,
                                                             find_sites_in_window(subsequence_coor,perm[0]),
                                                             find_sites_in_window(subsequence_coor,perm[1]))
                        
                        if not examin_each_comb_mass: #noneed eo examine each edited peptides mass - append peptide to peptides list
                            rna_pep_perm = rna_peptide(edited_sub_seq,subsequence_coor,a2g = perm[0], c2t = perm[1], cleavage_sites = perm_cs, edge = edge)
                            rna_peptides_list.append(rna_pep_perm)
#                            rna_pep_perm.print_peptide_data() 
                        
                        else: #each peptide mass should be examined before appended to list
#                            print(''.join(str(Seq(edited_sub_seq, generic_rna)).split('*')))
                            comb_mass = ProteinAnalysis(''.join(str(Seq(edited_sub_seq, generic_rna).translate()).split('*'))).molecular_weight()
#                            print(perm)
#                            print(comb_mass)
                            if comb_mass <= max_mass:
#                                print('in')
                                rna_pep_perm = rna_peptide(edited_sub_seq,subsequence_coor,a2g = perm[0], c2t = perm[1], cleavage_sites = perm_cs, edge = edge)
                                rna_peptides_list.append(rna_pep_perm)
#                                rna_pep_perm.print_peptide_data()
                            else:   #peptide mass is out of range - append to designated dictionary
                                if coor_str in over_mass_pep_combs_dict:
                                    over_mass_pep_combs_dict[coor_str].append(perm)
                                else:
                                    over_mass_pep_combs_dict.update({coor_str:[perm]})

# =============================================================================
            j+=1

    return rna_peptides_list, under_min_length, over_max_sites_peps, over_mass_peps, over_mass_pep_combs_dict


"""
given editing coordinations and a sequence and a lookahead lookbehind range
run over all codons, edit and if codos could be a cleavage site append to list.
create another binary list - for each codon fixed or mutable
"""
#a,b = create_all_cleavage_sites(sequence,digestion_rule,a2g,c2t,look_ahead_for_editnig,look_behind_for_editing)
def create_all_cleavage_sites(sequence,digestion_rule,a2g,c2t,look_ahead_for_editnig,look_behind_for_editing):
    
    cleavage_sites = deque([0])
    binary_for_fixed_sites = deque([1])
    length = len(sequence)
    
    #checking that sequence is devided to whol codons
    if len(sequence)%3:
        return None, None
    else:   
        for i in range(0,length-2,3): #run over all codons not including last wich is a fixed site
            
            #sequence is close to beginning, looking for permutations from beginning
            if i <= look_behind_for_editing:
                coor = (0,i+look_ahead_for_editnig-1)
                codon_to_asses_in_window = i
            #sequence is close to end, looking for permutations till end
            elif length-i<= look_ahead_for_editnig:
                coor = (i-look_behind_for_editing,length-1)
                codon_to_asses_in_window = look_behind_for_editing 
            else:
                coor = (i-look_behind_for_editing,i+look_ahead_for_editnig-1)
                codon_to_asses_in_window = look_behind_for_editing
            
            seq = (sequence[coor[0]:coor[1]+1])
            
            #all editing permutations in range
            permutations = find_permutations_in_range(creat_editing_sites_permutations_one_type(find_sites_in_window(coor,a2g)),
                                                      creat_editing_sites_permutations_one_type(find_sites_in_window(coor,c2t)))
            
            is_cs, is_fixed = determine_site(seq,coor,digestion_rule,codon_to_asses_in_window,permutations)
            
            if is_cs:
                cleavage_sites.append(i)
                binary_for_fixed_sites.append(is_fixed)
                
        cleavage_sites.append(length)
        binary_for_fixed_sites.append(1)
        
        return cleavage_sites,binary_for_fixed_sites
            
                
"""
for a subsequence determine codon_to_asses is a cleavage site and if so if mutable
"""        
#a,b = determine_site('AAACGC',[0,5],digestion_rule,3,permutations)
def determine_site(subsequence,coor,digestion_rule,codon_to_asses,permutations):
    
    cleavage_site = set()
    
    for perm in permutations:
        seq = edit_rna_as_peptide(subsequence,coor,perm[0],perm[1])
#        print(seq)
        cleavage_sites = find_cleavage_sites(digestion_rule,seq)
#        print(cleavage_sites)
        
        if codon_to_asses in cleavage_sites:
            cleavage_site.add(1)
        else:
            cleavage_site.add(0)
        
        if len(cleavage_site) == 2:
            break
        
    if len(cleavage_site)==2:
        return 1,0
    elif list(cleavage_site)[0] == 1:
        return 1,1
    else:
        return 0,1
       
            
"""
find all in-frame cleavage sites given a sequence and a digestion_rule
also append 0 and sequence_length if needed
"""
#k = find_cleavage_sites(digestion_rule,sequnece,append_zero=False, append_seq_end = False)
def find_cleavage_sites(digestion_rule,sequence,append_zero=False, append_seq_end = False):
    
    if append_zero:
        initial_list = [0] + [x.end() for x in digestion_rule.finditer(sequence)]
    else:
        initial_list = [x.end() for x in digestion_rule.finditer(sequence)]
    initial_set = set(initial_list)

    for site in initial_list:
        match1 = digestion_rule.match(sequence,site-1)
        match2 = digestion_rule.match(sequence,site-2)
        if match1:
            initial_set.add(match1.end())
        if match2:
            initial_set.add(match2.end()) 
 
    #return only sites representing in-frame codons
    cleavage_sites = [x for x in initial_set if not x%3 and x!=len(sequence)]
    
    if append_seq_end:
        cleavage_sites.append(len(sequence))
    
    cleavage_sites.sort()
    return cleavage_sites
            

"""
change a sub-sequence given:
the subsequence coordinates relative to the sequence
a_to_g editing coordinates relative to the sequence
c_to_t editing coordinates relative to the sequence
"""
#k = edit_rna_as_peptide(tsubseq, [273,314], [305,306], [])
def edit_rna_as_peptide(sequence, sequence_coordinates, a_to_g_list, c_to_t_list):
    splited_ag = [sequence[i:j] for i,j in zip([0] + sorted([x+1-sequence_coordinates[0] for x in a_to_g_list]), sorted([x-sequence_coordinates[0] for x in a_to_g_list]) + [None])]
    ag_swaped_seq = 'G'.join(splited_ag)
    splited_ct = [ag_swaped_seq[i:j] for i,j in zip([0] + sorted([x+1-sequence_coordinates[0] for x in c_to_t_list]), sorted([x-sequence_coordinates[0] for x in c_to_t_list]) + [None])]
    new_seq = 'T'.join(splited_ct)
    return new_seq

  

"""
return all combinations of elements from to lists
"""      
#k = find_permutations_in_range(creat_editing_sites_permutations_one_type(a2g_sites),creat_editing_sites_permutations_one_type(c2t_sites))
def find_permutations_in_range(a_to_g_permutations,c_to_t_permutations):
    return[(l[0],l[1]) for l in it.product(a_to_g_permutations,c_to_t_permutations)]
    
   
"""
return all permutation of a list
"""
def creat_editing_sites_permutations_one_type(lst):  
    editning_perm = []
    for l in range(len(lst)+1):
        [editning_perm.append(list(subset)) for subset in it.combinations(lst, l)]
    return editning_perm



"""
return a list of all elements from a given list that are within a given range
"""
#a = find_sites_in_window((18,56),[1,2,4,5,16,17,29,30,31,54])
def find_sites_in_window(editing_window,editing_sites):
    editing_sites_in_window = []
    editing_sites_it = it.chain(editing_sites)
    for i in editing_sites_it:
        if i<=editing_window[1]:
            if i >= editing_window[0]:
                editing_sites_in_window.append(i)
        else:
            break
    return editing_sites_in_window


def get_mass_range(sequence,coor,a2g_sites,c2t_sites):
    
    aa_seq_of_max_mass = []
    aa_seq_of_min_mass = []

    for i in range(0,len(sequence),3):
        
        max_aa_mass = 0
        min_aa_mass = 99999
        
        codon = sequence[i:i+3]
        editing_window = [i+coor[0], i+2+coor[0]]
        codon_a2g = find_sites_in_window(editing_window,a2g_sites)
        codon_c2t = find_sites_in_window(editing_window,c2t_sites)
        permutations = find_permutations_in_range(creat_editing_sites_permutations_one_type(codon_a2g),
                                                  creat_editing_sites_permutations_one_type(codon_c2t))
#        print('i=' + str(i))
        for perm in permutations:
            edited_sub_seq = edit_rna_as_peptide(codon,editing_window,
                                                 find_sites_in_window(editing_window,perm[0]),
                                                 find_sites_in_window(editing_window,perm[1]))
            
            aa = str(Seq(edited_sub_seq, generic_rna).translate())
            aa_mass = ProteinAnalysis(aa).molecular_weight()

            if aa_mass>max_aa_mass:
                max_aa_mass = aa_mass
                heavy_aa = aa
            if aa_mass<min_aa_mass:
                min_aa_mass = aa_mass
                light_aa = aa
                
        aa_seq_of_max_mass.append(heavy_aa)
        aa_seq_of_min_mass.append(light_aa)
        
    maximal_mass = ProteinAnalysis(''.join(aa_seq_of_max_mass)).molecular_weight()
    minimal_mass = ProteinAnalysis(''.join(aa_seq_of_min_mass)).molecular_weight()
    
    return minimal_mass, maximal_mass