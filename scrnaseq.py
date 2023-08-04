from data import *

class scRNAseq(Data):
    '''
    read scRNA-seq analysis performed by Seurat
    '''
    def __init__(self, config):
        self.config = config
        self.read_degs()
        self.read_vargenes()
        self.load_total_degs()
        self.read_clustern()
        #self.loadCellType()
        #self.loadCellRatio()
        
    def read_degs(self):
        self.CLUSTER_DEG = pd.read_csv(self.config['degs_path'], header=0, sep='\t')
        #self.CLUSTER_DEG['label'] = ['D' + str(i+1) for i in self.CLUSTER_DEG.cluster]
    
    def read_vargenes(self):
        self.VARGENE = pd.read_csv(self.config['vargene_path'], header=None, sep='\t')
        self.VARGENE.columns = ['gene', 'stdvar']
        
    def load_total_degs(self):
        '''
        number of scRNAseq cluster DEGs 
        '''
        overlap_list = []
        for i in range(np.max(self.CLUSTER_DEG.cluster)+1):
            len_overlap = self.CLUSTER_DEG[self.CLUSTER_DEG.cluster==i].shape[0]
            overlap_list.append(len_overlap)
        self.TOTAL_DEGS = overlap_list
    
    def read_clustern(self):
        self.CLUSTER_NUMBER = pd.read_csv(self.config['cluster_number_path'], header=None, sep='\t')
        self.CLUSTER_NUMBER.columns = ['cluster', 'number']
        self.CLUSTER_NUMBER['num_norm'] = self.CLUSTER_NUMBER.number / np.sum(self.CLUSTER_NUMBER.number)
    
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

