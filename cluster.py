import schicluster
from analysis import *
from genome import *

class Cluster(Analysis):
    '''
    Create clusters
    '''
    def __init__(self, config, genome='hg19', resolution=1000000):
        '''
        self.CELL: cell id, sample name and file path
        self.genome: an instance of Genome class
        self.SAMPLE_COUNT: list of sample type and their count
        '''
        ### init genome
        self.config = config
        self.genome = Genome(config = config, genome=genome, resolution=resolution)
        ### read cell path
        self.CELL = pd.read_csv(self.config['cell_path'], header=None, sep="\t")
        self.CELL.columns = ['cell_id', 'sample_type', 'path']
        self.CELL.sort_values(by=['sample_type', 'cell_id'], ascending=True, inplace=True)
        self.CELL.reset_index(drop = True, inplace = True)
        self.SAMPLE_COUNT = self.CELL.groupby('sample_type').size()
        
    def down_sampling(self, frac = 0.5):
        '''
        down sampling at a certain fraction, default fraction is 0.5.
        '''
        self.CELL = self.CELL.sample(frac = frac).copy()
        self.CELL.sort_values(by=['sample_type', 'cell_id'], ascending=True, inplace=True)
        self.CELL.reset_index(drop = True, inplace = True)
        self.SAMPLE_COUNT = self.CELL.groupby('sample_type').size()
        
    def read_contact(self):
        '''
        load contact matrix of single cell files and save to self.CONTACT_LIST
        '''
        file_name = self.CELL['path']
        self.CONTACT_LIST = []
        for i in range(len(file_name)):
            print('Read ' + str(i+1) + '/' + str(len(file_name)) + ' of single-cell contact files')
            contact_matrix = pd.read_csv(file_name[i], header=None, sep="\t")
            contact_matrix.columns = ["bin1", "bin2", "contact"]
            self.CONTACT_LIST.append(contact_matrix)
            
    def preprocessing(self, count_cutoff = 2):
        '''
        preprocess contact matrix for input of schicluster
        self.CELL_LIST: cells passed quality control
        count_cutoff must be > 1, otherwise hicluster_cpu will have an error
        '''
        self.CELL_LIST = []
        for i in range(len(self.CONTACT_LIST)):
            print('preprocessing ' + str(i+1) + '/' + str(len(self.CONTACT_LIST)))
            chr_counts = []
            for chr_number in self.genome.CHROMOSOME:
                max_num = max(self.genome.BINBED[self.genome.BINBED['chr'] == chr_number]['bin'])
                min_num = min(self.genome.BINBED[self.genome.BINBED['chr'] == chr_number]['bin'])
                chr_contact = self.CONTACT_LIST[i][(self.CONTACT_LIST[i]['bin1'] >= min_num) & (self.CONTACT_LIST[i]['bin1'] <= max_num) & (self.CONTACT_LIST[i]['bin2'] >= min_num) & (self.CONTACT_LIST[i]['bin2'] <= max_num)].copy()
                chr_contact = pd.concat([chr_contact.iloc[:,0:2] - min_num, chr_contact.iloc[:,2]], axis = 1)
                chr_counts.append(chr_contact.shape[0])
                chr_contact.to_csv(self.config['temp_path'] + '/' + self.CELL['cell_id'][i] + '_' + chr_number + '.txt', sep='\t', header=False, index=False)
            add_sample = True
            for j in chr_counts:
                if j < count_cutoff:
                    add_sample = False
            if add_sample:
                self.CELL_LIST.append(self.CELL['cell_id'][i])
    
    def call_pca(self):
        '''
        get self.RAW_PCA from contact matrix by schicluster
        '''
        network = self.CELL_LIST
        chromsize=self.genome.CHROM_SIZE
        nc=2
        res=self.genome.RESOLUTION
        ncpus=5
        os.chdir(self.config['temp_dir'])
        self.RAW_PCA = schicluster.hicluster_cpu(network=network, chromsize=chromsize, nc=nc, res=res, ncpus=ncpus)
        os.chdir(WORK_DIR)
    
    def save_pca(self):
        '''
        save self.RAW_PCA to output_dir as file raw_pca.txt
        '''
        os.chdir(self.config['output_dir'])
        np.savetxt('raw_pca_1.txt', self.RAW_PCA[0], fmt='%i')
        np.savetxt('raw_pca_2.txt', self.RAW_PCA[1])
        os.chdir(WORK_DIR)
        
    def read_pca(self):
        '''
        read self.RAW_PCA from output_dir
        '''
        os.chdir(self.config['output_dir'])
        raw_pca_1 = np.loadtxt('raw_pca_1.txt', dtype=int)
        raw_pca_2 = np.loadtxt('raw_pca_2.txt')
        self.RAW_PCA = (raw_pca_1, raw_pca_2)
        os.chdir(WORK_DIR)
        
    def call_cluster(self, X, n_cluster = 14):
        '''
        call cluster with number of cluster at n_cluster
        save cluster information to self.CLUSTER_LABEL
        save cluster information to file pdclsuter.txt
        '''
        scorelist = []
        for i in range(1, 1001, 1):
            kmeans = GaussianMixture(n_components=n_cluster,random_state=i).fit(X)
            score = silhouette_score(X, kmeans.predict(X))
            scorelist.append(score)
            print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))
        maxseed = np.argmax(pd.Series(scorelist))
        model = GaussianMixture(n_components=n_cluster, random_state=maxseed)
        model.fit(X)
        yhat = model.predict(X)
        self.CLUSTER_LABEL = yhat
        self.PDCLUSTER = pd.DataFrame({'cell_id': self.CELL.cell_id, 'sample_type': self.CELL.sample_type, 'cluster_id': yhat, 'PC1': X[:, 0], 'PC2': X[:, 1], 'PC3': X[:, 2], 'path': self.CELL.path})
        os.chdir(self.config['output_dir'])
        self.PDCLUSTER.to_csv('pdcluster.txt', index=False, header=True, sep='\t')
        print(self.PDCLUSTER.groupby(['cluster_id', 'sample_type']).size())
        os.chdir(WORK_DIR)
    
    def read_cluster(self):
        '''
        read pdcluster.txt from output_dir
        '''
        os.chdir(self.config['output_dir'])
        self.PDCLUSTER = pd.read_csv('pdcluster.txt', header=0, sep="\t")
        self.CLUSTER_LABEL = self.PDCLUSTER['cluster_id']
        print(self.PDCLUSTER.groupby(['cluster_id', 'sample_type']).size())
        os.chdir(WORK_DIR)
