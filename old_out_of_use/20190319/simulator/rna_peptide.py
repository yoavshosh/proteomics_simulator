import itertools as it

class rna_peptide:
    
    def __init__(self,sequence,coordinates, edited_sites, cleavage_sites = [], edge = None, C_terminus = 'no_change',N_terminus = 'no_change', cs_within_pep = 'no_change', sites_in_pep_range = 0, cancelled_cs = False):
        self.seq = sequence
        self.coo = coordinates
        self.edited_sites = edited_sites
        self.cleavage_sites = cleavage_sites
        self.edge = edge
        self.C_terminus = C_terminus
        self.N_terminus = N_terminus
        self.cs_within_pep = cs_within_pep
        self.sites_in_pep_range = sites_in_pep_range
        self.cancelled_cs_in_pep = cancelled_cs
    
    def print_peptide_data(self):
        print(self.seq + '| coordinates: ' + str(self.coo) + ' | edited_sites: ' + str(self.edited_sites)
                + ' | cleavage_sites (including boundaries): ' + str(self.cleavage_sites))
    
    
    