from HiC import *
from Genome import *
from TAD import *
from scRNAseq import *

class scHiC(HiC):
    def __init__(self, genome=Genome(), deg=scRNAseq()):
        '''
        self.allChrno: chr1, chr2, ... chrX
        self.Refseq: Genome Refseq bed file.
        self.getBinbed: generate self.binbed which is the list of genome loci for bins.
        '''
        #genome = Genome()
        #deg = scRNAseq()
        #pdb.set_trace()
        self.allChrno = genome.allchrno
        #pdb.set_trace()
        self.refReq = genome.refseq[genome.refseq.chrno.isin(self.allchrno)].copy()
        #pdb.set_trace()
        self.refReq.reset_index(drop = True, inplace = True)
        self.getBinbed()
        self.binbed = self.binbed[self.binbed.chr.isin(self.allchrno)].copy()
        self.binbed.reset_index(drop = True, inplace = True)
        self.deg = deg
