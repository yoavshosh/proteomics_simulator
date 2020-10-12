import itertools as it

class rna_peptide:
    
    def __init__(self,sequence,coordinates, a2g=[], c2t=[], cleavage_sites = [], edge = None, C_terminus = 'no_change',N_terminus = 'no_change', cs_within_pep = 'no_change', sites_in_pep_range = 0, cancelled_cs = False):
        self.seq = sequence
        self.coo = coordinates
        self.a2g = a2g
        self.c2t = c2t
        self.cleavage_sites = cleavage_sites
        self.edge = edge
        self.C_terminus = C_terminus
        self.N_terminus = N_terminus
        self.cs_within_pep = cs_within_pep
        self.sites_in_pep_range = sites_in_pep_range
        self.cancelled_cs_in_pep = cancelled_cs
    
    def print_peptide_data(self):
        print(self.seq + '| coordinates: ' + str(self.coo) + ' | edited_sites: ' + str(self.a2g)+';'+str(self.c2t) 
                + ' | cleavage_sites (including boundaries): ' + str(self.cleavage_sites))
    

    def update_a2g_editing_sites(self, editing_sites):
        self.a2g = []
        editing_sites_it = it.chain(editing_sites)
        for i in editing_sites_it:
            if i<=self.coo[1]:
                if i >= self.coo[0]:
                    self.a2g.append(i)
            else:
                break

    def update_c2t_editing_sites(self, editing_sites):
        self.c2t = []
        editing_sites_it = it.chain(editing_sites)
        for i in editing_sites_it:
            if i<=self.coo[1]:
                if i >= self.coo[0]:
                    self.c2t.append(i)
            else:
                break

    
    def find_cleavage_sites(self, digestion_rule, cleave_before_pattern=False,append_zero = False, append_seq_end = False):
                    
        if append_zero: #peptide is at the beginning of sequence thus appending a pseudo cleavage site - 0
            initial_list = [0] + [x.end() for x in digestion_rule.finditer(self.seq)]
        else:
            initial_list = [x.end() for x in digestion_rule.finditer(self.seq)]
        initial_set = set(initial_list)
    
        #find other overlapping sites
        if cleave_before_pattern: #cleavage sites are before regex pattern (usually not the casee)
            for site in initial_list:
                match1 = digestion_rule.match(self.seq,site+1)
                match2 = digestion_rule.match(self.seq,site+2)
                if match1:
                    initial_set.add(match1.start())
                if match2:
                    initial_set.add(match2.start())
                
        else:
            for site in initial_list:
                match1 = digestion_rule.match(self.seq,site-1)
                match2 = digestion_rule.match(self.seq,site-2)
                if match1:
                    initial_set.add(match1.end())
                if match2:
                    initial_set.add(match2.end()) 
        
        #return only sites representing in-frame codons
        cleavage_sites = [x + self.coo[0] for x in initial_set if not x%3 and x!=len(self.seq)]
        
        if append_seq_end: #peptide is at the end of sequence thus appending a pseudo cleavage site - 
            cleavage_sites.append(len(self.seq)+self.coo[0])
        cleavage_sites.sort()
    
        self.cleavage_sites = cleavage_sites
        
    
         
        
        
    