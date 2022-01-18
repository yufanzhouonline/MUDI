from HiC import *
from Genome import *
from TAD import *
from scRNAseq import *


class scHiC(HiC):
    def __init__(self, genome=Genome(), deg=scRNAseq()):
        '''
        self.allChrno: chr1, chr2, ... chrX
        self.refSeq: Genome Refseq bed file.
        self.getBinBed: generate self.binbed which is the list of genome loci for bins.
        '''
        #genome = Genome()
        #deg = scRNAseq()
        #pdb.set_trace()
        self.allChrNo = genome.allChrNo
        #pdb.set_trace()
        self.refSeq = genome.refSeq[genome.refSeq.chrNo.isin(self.allChrNo)].copy()
        #pdb.set_trace()
        self.refSeq.reset_index(drop = True, inplace = True)
        self.getBinBed()
        self.binBed = self.binBed[self.binBed.chr.isin(self.allChrNo)].copy()
        self.binBed.reset_index(drop = True, inplace = True)
        self.deg = deg
        
