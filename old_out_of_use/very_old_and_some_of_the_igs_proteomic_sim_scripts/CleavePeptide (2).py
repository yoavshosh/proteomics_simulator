import re
import itertools as it
from collections import deque


"""
    Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.

        .. note::
            The sequence is expected to be in one-letter uppercase notation.
            Otherwise, some of the cleavage rules in :py:data:`expasy_rules`
            will not work as expected.

    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
        :py:data:`expasy_rules` contains cleavage rules for popular cleavage agents.
    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.

        ..note ::
            This checks for string length, which is only correct for one-letter
            notation and not for full *modX*. Use :py:func:`length` manually if
            you know what you are doing and apply :py:func:`cleave` to *modX*
            sequences.

    Returns
    -------
    out : set
        A set of unique (!) peptides.
        
"""

def cleave(sequence, rule, missed_cleavages=0, min_length=None):

    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
    for i in it.chain([x.end() for x in re.finditer(rule, sequence)],
                      [None]):
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append(seq)
    return peptides



#==============================================================================
# """
# this function cleaves sequences of amino acids
# input parameters:
# sequence : str - The sequence of a polypeptide.
# rule : str - A string with a regular expression describing the C-terminal site of cleavage.    
# missed_cleavages : int, optional - The maximal number of allowed missed cleavages. Defaults to 0.
# overlap : bool, optional - 
#         Set this to :py:const:`True` if the cleavage rule is complex and
#         it is important to get all possible peptides when the matching
#         subsequences overlap (e.g. 'XX' produces overlapping matches when
#         the sequence contains 'XXX'). Default is :py:const:`False`.
#         Use with caution: enabling this results in exponentially growing
#         execution time.
# 
# output: A set of unique (!) peptides.
# Note - this function was taken from Pyteomics by Lev Levitsky 
# """
# 
# 
# def cleave(sequence, rule, missed_cleavages=0, overlap=False):
#     
#     peptides = set()
#     cleavage_sites = deque([0], maxlen=missed_cleavages+2)
#     
#     for i in chain(map(lambda x: x.end(), re.finditer(rule, sequence)),[None]):
#         cleavage_sites.append(i)
#         for j in range(0, len(cleavage_sites)-1):
#             peptides.add(sequence[cleavage_sites[j]:cleavage_sites[-1]])
#         if overlap and i not in {0, None}:
#             peptides.update(cleave(sequence[i:], rule, missed_cleavages, overlap))
# 
#     if '' in peptides:
#         peptides.remove('')
#         
#     return peptides
#==============================================================================


"""
expasy_rules = {
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
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
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
