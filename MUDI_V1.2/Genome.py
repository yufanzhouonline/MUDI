from NGS import *

class Genome(object):
    '''
    Class of Genome.
    '''
    def __init__(self, genome='hg19'):
        if genome=='hg19':
            self.chrNo()
            self.loadRefSeq()
            self.geneClass()
        else:
            print('ERROR! Please input the correct genome.')
    
    def chrNo(self):
        self.allChrNo = ["chr" + str(i+1) for i in range(22)] + ['chrX']
    
    def loadRefSeq(self, refSeqPath = './data/hg19/refseq.bed'):
        refSeq = pd.read_csv(refSeqPath, header=None, sep='\t')
        #refSeq.columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
        #    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score',
        #    'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
        refSeq.columns = ['chrom', 'txStart', 'txEnd', 'name2', 'exon', 'strand']
        self.refSeq = refSeq.loc[:,['name2', 'chrom', 'strand', 'txStart', 'txEnd']].copy()
        self.refSeq.drop_duplicates(subset=None, keep='first', inplace=True)
        self.refSeq.reset_index(drop = True, inplace = True)
        self.refSeq.columns = ['gene', 'chrNo', 'strand', 'start', 'end']
    
    def geneClass(self, g1sPath = './data/signaling/g1s.txt', g2mPath = './data/signaling/g2m.txt', cellCyclePath = './data/signaling/cellcycle.txt', houseKeepingPath = './data/signaling/housekeeping.txt'):
        self.g1s = pd.read_csv(g1sPath, header=None, sep="\t")
        self.g1s.columns = ["gene"]
        self.g2m = pd.read_csv(g2mPath, header=None, sep="\t")
        self.g2m.columns = ["gene"]
        self.cellcycle = pd.read_csv(cellCyclePath, header=None, sep="\t")
        self.cellcycle.columns = ["gene"]
        self.housekeeping = pd.read_csv(houseKeepingPath, header=None, sep="\t")
        self.housekeeping.columns = ["gene"]
