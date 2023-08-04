from data import *

class Genome(Data):
    '''
    Class of Genome.
    '''
    def __init__(self, config, genome='hg19', resolution=1000000):
        ### init resolution
        self.config = config
        self.RESOLUTION = resolution
        ### init genome
        if genome == 'hg19':
            self.chr_number()
            self.read_refseq()
            self.gene_class()
            self.read_binbed()
        else:
            print('ERROR! Please input the correct genome.')
    
    def chr_number(self, genome='hg19'):
        if genome == 'hg19':
            #self.CHROMOSOME = ["chr" + str(i+1) for i in range(22)] + ['chrX'] + ['chrY'] + ['chrM']
            self.CHROMOSOME = ["chr" + str(i+1) for i in range(22)] + ['chrX']
            hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
            chrom = [str(i+1) for i in range(22)] + ['X']
            self.CHROM_SIZE = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
    
    def read_refseq(self):
        REFSEQ = pd.read_csv(self.config['refseq_path'], header=None, sep='\t')
        #refSeq.columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
        #    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score',
        #    'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
        REFSEQ.columns = ['chrom', 'txStart', 'txEnd', 'name2', 'exon', 'strand']
        self.REFSEQ = REFSEQ.loc[:,['name2', 'chrom', 'strand', 'txStart', 'txEnd']].copy()
        self.REFSEQ.drop_duplicates(subset=None, keep='first', inplace=True)
        self.REFSEQ = self.REFSEQ[self.REFSEQ.chrom.isin(self.CHROMOSOME)].copy()
        self.REFSEQ.reset_index(drop = True, inplace = True)
        self.REFSEQ.columns = ['gene', 'chr', 'strand', 'start', 'end']
    
    def read_binbed(self):
        '''
        self.BINBED: list of genome loci for bins in a certain resolution.
        '''
        if self.RESOLUTION == 1000000:
            bed_path = self.config['bed_path']
        self.BINBED = pd.read_csv(bed_path, header=None, sep='\t')
        self.BINBED.columns = ['chr', 'start', 'end', 'bin']
        self.BINBED = self.BINBED[self.BINBED.chr.isin(self.CHROMOSOME)].copy()
    
    def gene_class(self):
        '''
        define gene class: G1S, G2M, CELLCYCLE and HOUSEKEEPING genes
        '''
        self.G1S = pd.read_csv(self.config['g1s_path'], header=None, sep="\t")
        self.G1S.columns = ["gene"]
        self.G2M = pd.read_csv(self.config['g2m_path'], header=None, sep="\t")
        self.G2M.columns = ["gene"]
        self.CELLCYCLE = pd.read_csv(self.config['cellcycle_path'], header=None, sep="\t")
        self.CELLCYCLE.columns = ["gene"]
        self.HOUSEKEEPING = pd.read_csv(self.config['housekeeping_path'], header=None, sep="\t")
        self.HOUSEKEEPING.columns = ["gene"]
