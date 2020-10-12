
#path of souces file AS A RAW STRING!!!! (the letter r before the string, exmp: r'the_path')
father_path = r'C:\Users\user\Google_Drive\NGS_project\Proteomics_project\files\source_sequences\from_yariv'

#name of sources file (including '.fasta' ending)
input_file = r'simulation.fasta'

# a string in records description before clone id
clone_identifier = r'CDR3: '

#True if CDR3 elements will include also C-terminus untill sequnces end
continue_CDR3_element = False

#lower and upper boound of molecular weight threshold as [a,b] where a is lower bound.
#[] (empty list) is now bounds -> no molecular weight filtration for peptides
mw_threshold_filtration = [350,6000] 


#list of numbers of iss cleavages allowed in peptides.
#for exeple as [1,2] filter only peptides that contain 1 or 2 miss cleavages.
#[] (empty list) -> no miss-cleavages filtration for peptides
miss_cleavages_filtration = []

# True if filtration process should also yield a fasta file containing only the filtered peptides
# Note that all of the peptides are present in the unfiltered databases  anyway
# and molecular weight and missed cleavage sites could be easaly deduced using the function "python peptides_analysis.py peptide_data".
# choosing False may save running time and storage.
create_filtered_database = True

#This are some cleavage rules in regular expressions
#use may use this dictionary to call the rule for each iteration
rules = {'trypsin':r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
         'asp-n':r'\w(?=D)',
         'lys-c':r'K(?!$)',
         'arg-c':r'R(?!$)',
         'glutamyl endopeptidase': r'E(?!$)',
         'chymotrypsin high specificity':r'([FY](?=[^P]))|(W(?=[^MP]))'}

"""
expacy_rules = {
    'arg-c':         r'R',
    'asp-n':         r'\w(?=D)',
    'bnps-skatole' : r'W',
    'caspase 1':     r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2':     r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3':     r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4':     r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5':     r'(?<=[LW]EH)D',
    'caspase 6':     r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7':     r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8':     r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9':     r'(?<=LEH)D',
    'caspase 10':    r'(?<=IEA)D',
    'chymotrypsin high specificity' : r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity': r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   r'R',
    'cnbr':          r'M',
    'enterokinase':  r'(?<=[DE]{3})K',
    'factor xa':     r'(?<=[AFGILTVM][DE]G)R',
    'formic acid':   r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b':    r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc':          r'K',
    'ntcb':          r'\w(?=C)',
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k':  r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin':   r'[^DE](?=[AFILMV])',
    'thrombin':      r'((?<=G)R(?=G))|'
                     r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
    }
    
"""


#ruls list - a list of "sublists" of ruls.
#each sublists represents proteomics using certain enzyme(s) with pre-determined efficiency(ies)
#each list in the list contains  number of dictionaries
#each dictionary represents a cleavage rule and efficiency
#for example: [{'name':'lys-c','rule':rules['lysc'],'miss_cleavages':0,'min_length':None}]
#means that ther is only one enzyme (one dictionary in the list)
#name of rule is lys-c. the regex rule is taken from the rules dictionary above, max miss cleavages allowed for peptide = 0, no min length for peptides
#if more than 1 dictionaty in a list -> proteomics with more than one enzyme
rules_list = [
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':1,'min_length':None}],
              [{'name':'trypsin','rule':rules['trypsin'],'miss_cleavages':1,'min_length':None}],
              [{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':1,'min_length':None}],
              [{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':1,'min_length':None}],
              [{'name':'arg-c','rule':rules['arg-c'],'miss_cleavages':1,'min_length':None}],
              [{'name':'chymotrypsin','rule':rules['chymotrypsin high specificity'],'miss_cleavages':1,'min_length':None}],
              [{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':1,'min_length':None},{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':1,'min_length':None}],
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':1,'min_length':None},{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':1,'min_length':None}],
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':1,'min_length':None},{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':1,'min_length':None}],
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':0,'min_length':None}],
              [{'name':'trypsin','rule':rules['trypsin'],'miss_cleavages':0,'min_length':None}],
              [{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':0,'min_length':None}],
              [{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':0,'min_length':None}],
              [{'name':'arg-c','rule':rules['arg-c'],'miss_cleavages':0,'min_length':None}],
              [{'name':'chymotrypsin','rule':rules['chymotrypsin high specificity'],'miss_cleavages':0,'min_length':None}],
              [{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':0,'min_length':None},{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':0,'min_length':None}],
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':0,'min_length':None},{'name':'lys-c','rule':rules['lys-c'],'miss_cleavages':0,'min_length':None}],
              [{'name':'asp-n','rule':rules['asp-n'],'miss_cleavages':0,'min_length':None},{'name':'glu-c','rule':rules['glutamyl endopeptidase'],'miss_cleavages':0,'min_length':None}]
             ]




