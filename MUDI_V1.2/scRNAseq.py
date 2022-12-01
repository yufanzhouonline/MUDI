from RNAseq import *

class scRNAseq(RNAseq):
    def __init__(self):
        self.loadClusterDEG()
        self.loadVarGene()
        self.loadTotalDEG()
        self.loadClusterN()
        #self.loadCellType()
        #self.loadCellRatio()
        
    def loadClusterDEG(self, clusterDEGpath = './data/samples/clusterDEGsnew.txt'):
        self.clusterDEG = pd.read_csv(clusterDEGpath, header=0, sep='\t')
        self.clusterDEG['label'] = ['D' + str(i+1) for i in self.clusterDEG.cluster]
    
    def loadVarGene(self, varGenePath = './data/samples/varGeneNew.rnk'):
        self.varGene = pd.read_csv(varGenePath, header=None, sep='\t')
        self.varGene.columns = ['gene', 'stdvar']
        
    def loadTotalDEG(self):
        '''
        number of scRNAseq cluster DEGs 
        '''
        overlapList = []
        for i in range(np.max(self.clusterDEG.cluster)+1):
            lenOverlap = self.clusterDEG[self.clusterDEG.cluster==i].shape[0]
            overlapList.append(lenOverlap)
        self.totalDEG = overlapList
    
    def loadClusterN(self, clusterNumberPath = './data/samples/clusterNumber.txt'):
        self.clusterNumber = pd.read_csv(clusterNumberPath, header=None, sep='\t')
        self.clusterNumber.columns = ['cluster', 'number']
        self.clusterNumber['numNorm'] = self.clusterNumber.number / np.sum(self.clusterNumber.number)
    
    '''
    def loadCellType(self, cellTypePath = './data/samples/clusterCellType.txt'):
        ###cell type of scRNA-seq clusters
        self.cellType = pd.read_csv(cellTypePath, header=0, sep='\t')
        self.cellType.columns = ['label', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12', 'D13']
        self.cellType['cell'] = ['MCF7', 'MCF7', 'MCF7M1', 'MCF7M1', 'MCF7TR', 'MCF7TR']
    
    def loadCellRatio(self):
        self.cellRatio = pd.DataFrame(self.cellType[self.cellType.cell == 'MCF7'].sum(axis=0).iloc[1:14])
        self.cellRatio.columns = ['mcf7']
        self.cellRatio['mcf7m1']= self.cellType[self.cellType.cell == 'MCF7M1'].sum(axis=0).iloc[1:14]
        self.cellRatio['mcf7tr']= self.cellType[self.cellType.cell == 'MCF7TR'].sum(axis=0).iloc[1:14]
        self.cellRatio['maxcell'] = self.cellRatio.max(axis=1)
        self.cellRatio['total'] = pd.DataFrame(self.cellType.sum(axis=0).iloc[1:14])
        self.cellRatio['ratio'] = self.cellRatio['maxcell'] / self.cellRatio['total']
        self.cellRatio.index = range(13)
    '''

