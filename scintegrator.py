from analysis import *
from genome import *
from scrnaseq import *

class Integrator(Analysis):
    def __init__(self, config):
        '''
        self.GENOME: instance of class Genome
        self.DEG: instance of class scRNAseq
        self.CHROMOSOME: chr1, chr2, ... chrX
        self.REFSEQ: Genome REFSEQ bed file.
        self.read_binbed: generate self.BINBED which is the list of genome loci for bins.
        '''
        self.config = config
        self.GENOME = Genome(config = config)
        self.DEG = scRNAseq(config = config)
        self.CHROMOSOME = self.GENOME.CHROMOSOME
        self.REFSEQ = self.GENOME.REFSEQ.copy()
        self.REFSEQ.reset_index(drop = True, inplace = True)
        self.read_binbed()
        
    def read_binbed(self):
        '''
        Read bin bed file.
        self.BINBED: list of genome loci for bins.
        '''
        self.BINBED = pd.read_csv(self.config['bed_path'], header=None, sep='\t')
        self.BINBED.columns = ['chr', 'bin1', 'bin2', 'pos']
        self.BINBED = self.BINBED[self.BINBED.chr.isin(self.CHROMOSOME)].copy()
        self.BINBED.reset_index(drop = True, inplace = True)
        
    def cell_label(self):
        '''
        Load the file of cell labeling.
        self.CELL: list of sample
        '''
        self.CELL = pd.read_csv(self.config['label_path'], header=None, sep="\t")
        self.CELL.columns = ['cell_id', 'sample_type', 'path']
    
    def sample_count(self, sample_list=[], label_list=[]):
        '''
        Count the samples and make the plot when there are three sample types.
        '''
        if not sample_list:
            sample_list = list(set(self.CELL['sample_type']))
        sample_counts = []
        sample_counts.append(self.CELL[self.CELL['sample_type'] == sample_list[0]].shape[0])
        sample_counts.append(self.CELL[self.CELL['sample_type'] == sample_list[1]].shape[0])
        sample_counts.append(self.CELL[self.CELL['sample_type'] == sample_list[2]].shape[0])
        ###Plot the sample count of single cell Hi-C
        if not label_list:
            label_list = sample_list
        rects1 = plt.bar(label_list, sample_counts, color = ['#88224488', '#22884488', '#22448888'])
        
        def autolabel(rects):
            for rect in rects:
                height = rect.get_height()
                #plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, '%s' % float(height))
                plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, str(height))
        
        autolabel(rects1)
        plt.ylim(0, 125)
        plt.ylabel('# of cells', fontsize=15)
        plt.show()
    
    def cell_cluster(self, label_prefix = 'C'):
        '''
        Load cell cluster file.
        self.PDCLUSTER: list of clusters.
        self.CLUSTER_LIST: list of sample name in each cluster
        self.CLUSTER_NUM: number of cluster type
        '''
        self.PDCLUSTER = pd.read_csv(self.config['cluster_path'], sep='\t', header=0)
        self.PDCLUSTER.columns = ['sample', 'cell', 'cluster', 'PC1', 'PC2', 'PC3', 'path']
        self.PDCLUSTER['label'] = [label_prefix + str(i+1) for i in self.PDCLUSTER.cluster]
        self.CLUSTER_LIST = []
        self.CLUSTER_NUM = len(set(self.PDCLUSTER.cluster))
        for i in range(self.CLUSTER_NUM):
            self.CLUSTER_LIST.append(list(self.PDCLUSTER[self.PDCLUSTER.cluster==i]['sample']))
    
    def common_binnum(self, sample_name, cluster_list, binvalue_cutoff = 0):
        '''
        Return the number of Commonly Associating Domain (CADs) of the cluster.
        '''
        group_name = sample_name[sample_name.cell_id.isin(cluster_list)].copy()
        group_name.reset_index(drop = True, inplace = True)
        set_list = []
        #res = 2000000
        for i in group_name.index:
            run_data = group_name.loc[i, 'path']
            print("Reading " + run_data + "......")
            csv_mat = pd.read_csv(run_data, header=None, sep="\t")
            csv_mat.columns = ["bin1", "bin2", "binvalue"]
            csv_mat = csv_mat[csv_mat.binvalue > binvalue_cutoff]
            csv_mat['twobin'] = csv_mat.bin1 * 10000+ csv_mat.bin2
            set_list.append(set(csv_mat.twobin))
        common_set = set_list[0]          ###get the common_set of all cells to common_set
        for i in range(1,len(set_list)):
            common_set = common_set.intersection(set_list[i])
        ###calculate the common set
        common_set = sorted(list(common_set))
        ###add the bins of chr22:2846-2897
        #common_set.extend([i*10000+i for i in range(2846,2898)])
        #print(len(common_set))
        ###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
        ###
        return len(common_set)
    
    def common_bin_per_cell(self):
        '''
        Return CADs number per cell of all clusters.
        '''
        per_cell = []
        for i in range(self.CLUSTER_NUM):
            cell_number = self.PDCLUSTER[self.PDCLUSTER.cluster==i].shape[0]
            per_cell.append(self.common_binnum(self.CELL, self.CLUSTER_LIST[i], 1000000) / cell_number)
        return per_cell
    
    def common_binset(self, sample_name, cluster_list, binvalue_cutoff = 0):
        '''
        The set of CADs.
        '''
        group_name = sample_name[sample_name.cell_id.isin(cluster_list)].copy()
        group_name.reset_index(drop = True, inplace = True)
        set_list = []
        #res = 2000000
        for i in group_name.index:
            run_data = group_name.loc[i, 'path']
            print("Reading " + run_data + "......")
            csv_mat = pd.read_csv(run_data, header=None, sep="\t")
            csv_mat.columns = ["bin1", "bin2", "binvalue"]
            csv_mat = csv_mat[csv_mat.binvalue > binvalue_cutoff]
            csv_mat['twobin'] = csv_mat.bin1 * 10000+ csv_mat.bin2
            set_list.append(set(csv_mat.twobin))
        common_set = set_list[0]                ###get the common_set of all cells to common_set
        for i in range(1,len(set_list)):
            common_set = common_set.intersection(set_list[i])
        ###calculate the common set
        common_set = sorted(list(common_set))
        ###add the bins of chr22:2846-2897
        #common_set.extend([i*10000+i for i in range(2846,2898)])
        #print(len(common_set))
        ###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
        ###
        return common_set
    
    def cads_score(self, sample_name, cluster_list, binvalue_cutoff = 0):
        '''
        Calculate the score of CADs.
        self.BIN_LIST: list of interaction frequencies of each cell in a cluster.
        self.BIN_SCORE: data frame of interaction frequencies, common bins (# of cells) of a cluster.
        '''
        group_name = sample_name[sample_name.cell_id.isin(cluster_list)].copy()
        group_name.reset_index(drop = True, inplace = True)
        ###Prepare to get the list of clusters and store in the list of self.BIN_LIST
        self.BIN_LIST = []
        for i in group_name.index:
            run_data = group_name.loc[i, 'path']
            print("Reading " + run_data + "......")
            csv_mat = pd.read_csv(run_data, header=None, sep="\t")
            csv_mat.columns = ["bin1", "bin2", "binvalue"]
            csv_mat = csv_mat[csv_mat.binvalue > binvalue_cutoff]
            csv_mat['twobin'] = csv_mat.bin1 * pow(10, len(str(np.max(csv_mat.bin2)))) + csv_mat.bin2
            self.BIN_LIST.append(csv_mat)
        ###Generate the cluster bin with the number of cells in column of cell
        cluster_bin = self.BIN_LIST[0]
        cluster_bin['cell'] = 1
        for i in range(len(self.BIN_LIST))[1:]:
            print(i+1, '/', len(self.BIN_LIST))
            nextbin = self.BIN_LIST[i]
            nextbin['cell'] = 1
            common_list = list(set(cluster_bin.twobin).intersection(set(nextbin.twobin)))
            firstdiff = cluster_bin[~cluster_bin.twobin.isin(common_list)].copy()
            nextdiff = nextbin[~nextbin.twobin.isin(common_list)].copy()
            ###reset the cell number of both common bin
            both_common = cluster_bin[cluster_bin.twobin.isin(common_list)].copy()
            both_common['cell'] = both_common['cell'] + 1
            ###reset the bin value
            both_common2 = nextbin[nextbin.twobin.isin(common_list)].copy()
            both_common.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            both_common2.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            both_common['binvalue2'] = list(both_common2['binvalue'])
            both_common['binvalue'] = both_common.apply(lambda row : max(row['binvalue'], row['binvalue2']), axis=1)
            both_common.drop('binvalue2', axis=1, inplace=True)
            ###combine three df
            cluster_bin = pd.concat([firstdiff, nextdiff, both_common], axis = 0)
            cluster_bin.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            cluster_bin.reset_index(drop = True, inplace = True)
        ###save the score to binscore
        self.BIN_SCORE = cluster_bin.copy()
    
    def cads_loci(self, bed_df, score_df):
        '''
        Return the genes of CADs.
        beddf: self.BINBED
        score_df: self.BIN_SCORE
        '''
        ###get bins within binbed
        bin_list = list(self.BINBED.pos)
        score_df = score_df[score_df.bin1.isin(bin_list) & score_df.bin2.isin(bin_list)].copy()
        score_df.reset_index(drop = True, inplace = True)
        ###get the genome loci
        score_df['bin1chr'] = ''
        score_df['bin1start'] = 0
        score_df['bin1end'] = 0
        score_df['bin2chr'] = ''
        score_df['bin2start'] = 0
        score_df['bin2end'] = 0
        for i in score_df.index:
            print(i+1, '/', score_df.shape[0])
            score_df.loc[i, 'bin1chr']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin1']]['chr'])[0]
            score_df.loc[i, 'bin1start']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin1']]['bin1'])[0]
            score_df.loc[i, 'bin1end']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin1']]['bin2'])[0]
            score_df.loc[i, 'bin2chr']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin2']]['chr'])[0]
            score_df.loc[i, 'bin2start']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin2']]['bin1'])[0]
            score_df.loc[i, 'bin2end']= list(bed_df[bed_df.pos==score_df.loc[i, 'bin2']]['bin2'])[0]
        return score_df
    
    def cluster_gene(self, expanded_bin, mode='max'):
        '''
        Make genes frequencies of clusters from expanded bins which have genome loci.
        '''
        gene_cell = self.REFSEQ.copy()
        ###initiate cell number and interaction frequency
        gene_cell['cell'] = 0
        gene_cell['binvalue'] = 0
        for i in gene_cell.index:
            ### [a.start,a.end] and [b.start,b.end] overlap?
            ### if(A.start <= B.end && B.start <= A.end) overlap
            ### if(A.start > B.end || B.start > A.end)   No overlap
            gene_bin = expanded_bin[((gene_cell.loc[i, 'chr'] == expanded_bin.bin1chr) & (gene_cell.loc[i, 'start'] <= expanded_bin.bin1end) & (expanded_bin.bin1start <= gene_cell.loc[i, 'end'])) |
                                  ((gene_cell.loc[i, 'chr'] == expanded_bin.bin2chr) & (gene_cell.loc[i, 'start'] <= expanded_bin.bin2end) & (expanded_bin.bin2start <= gene_cell.loc[i, 'end']))].copy()
            print(i+1, '/', gene_cell.shape[0], gene_cell.loc[i, 'gene'])
            #print(genebin)
            if gene_bin.shape[0]!=0:
                if mode=='max':
                    gene_cell.loc[i, 'cell'] = np.max(gene_bin.cell)
                    gene_cell.loc[i, 'binvalue'] = np.max(gene_bin.binvalue)
                elif mode=='mean':
                    gene_cell.loc[i, 'cell'] = np.mean(gene_bin.cell)
                    gene_cell.loc[i, 'binvalue'] = np.mean(gene_bin.binvalue)
                elif mode=='min':
                    gene_cell.loc[i, 'cell'] = np.min(gene_bin.cell)
                    gene_cell.loc[i, 'binvalue'] = np.min(gene_bin.binvalue)
                else:
                    print("Wrong mode.")
                ###Normalize the binvalue
                gene_cell.loc[i, 'binvalue'] = np.log2(gene_cell.loc[i, 'binvalue'])
        gene_cell.sort_values(by=['gene', 'cell'], ascending=False, inplace=True)
        gene_cell = gene_cell.drop_duplicates(subset='gene', keep='first', inplace=False).copy()
        gene_cell.sort_values(by='gene', ascending=True, inplace=True)
        gene_cell.reset_index(drop = True, inplace = True)
        return gene_cell
    
    def cluster_gene_list(self, cutoff_list = []):
        '''
        Return the list of genes in CADs of each cluster to:
        self.GENE_LIST
        '''
        # reset cutoff_list as [0, 0, 0, ...]
        if not cutoff_list:
            cutoff_list = [0] * len(self.CLUSTER_LIST)
        ### get list of genes in CADs
        self.GENE_LIST = []
        for i in range(len(self.CLUSTER_LIST)):
            self.cads_score(sample_name = self.CELL, cluster_list = self.CLUSTER_LIST[i], binvalue_cutoff = cutoff_list[i])
            cads_bin = self.BIN_SCORE[self.BIN_SCORE.cell==len(self.CLUSTER_LIST[i])].copy()
            cads_bin.reset_index(drop = True, inplace = True)
            cads_bin = self.cads_loci(self.BINBED, cads_bin)
            gene_cell = self.cluster_gene(cads_bin, mode='max')
            self.GENE_LIST.append(gene_cell[gene_cell.cell==len(self.BIN_LIST)].copy())
    
    def integration_score(self):
        '''
        Return the gene list with integration score of each cluster
        allPdList[i][j]: i: cluster of scHi-C, j: cluster of scRNA-seq
        '''
        allpd_list = []
        for i in range(len(self.CLUSTER_LIST)):
            cluster_gene_cell = self.GENE_LIST[i]
            genepd_list = []
            for j in range(np.max(self.DEG.CLUSTER_DEG.cluster)+1):
                #pdb.set_trace()
                cluster_deg2 = self.DEG.CLUSTER_DEG[self.DEG.CLUSTER_DEG.cluster==j].copy()
                overlapped_gene = set(cluster_gene_cell.gene).intersection(set(cluster_deg2.gene))
                genelist_sorted = sorted(list(overlapped_gene))
                ###integration score
                listfg = cluster_gene_cell[cluster_gene_cell.gene.isin(overlapped_gene)].sort_values(by='gene', ascending=True, inplace=False).binvalue.copy()
                listeg = cluster_deg2[cluster_deg2.gene.isin(overlapped_gene)].sort_values(by='gene', ascending=True, inplace=False).avg_log2FC.copy()
                listfg.reset_index(drop = True, inplace = True)
                listeg.reset_index(drop = True, inplace = True)
                integration_gene = listfg * listeg / self.DEG.TOTAL_DEGS[j] / self.DEG.CLUSTER_NUMBER.num_norm[j] * np.sum(self.DEG.TOTAL_DEGS)
                if str(integration_gene)=='nan':
                    integration_gene = 0
                genepd = pd.DataFrame({'gene': genelist_sorted, 'integration': integration_gene})
                genepd.sort_values(by='integration', ascending=False, inplace=True)
                genepd_list.append(genepd)
            allpd_list.append(genepd_list)
        return allpd_list
    
    def common_binbed(self, sample_name, cluster_list):
        '''
        Return the bed format of CADs.
        chr: chromosome, bin1: start, bin2: end, pos: bin number.
        self.BINBED: bed file to parse genome loci of the bin.
        '''
        common_set = self.common_binset(sample_name = sample_name, cluster_list = cluster_list)
        ###loci of CADs
        loci1 = pd.Series(common_set)//10000
        loci2 = pd.Series(common_set)%10000
        for i in loci1.index:
            print(i)
            if loci1[i] == loci2[i]:
                print('True')
            else:
                print('False')
        common_bin_bed = self.BINBED[self.BINBED.pos.isin(loci1)].copy()
        return common_bin_bed
    
    def cell_shift(self, sample_name, pdcluster, clusterno, common_bin, res=100000, cads=True):
        '''
        Get the shifted TADs within CADs.
        clusterno: Number of cluster
        common_bin: Common bin of cluster (CADs)
        res: Resolution for TADs
        cads: True, CADs, False: NADs 
        '''
        load_sample = sample_name[sample_name.cell_id.isin(pdcluster[pdcluster.cluster==clusterNo]['sample'])].copy()
        load_sample.reset_index(drop = True, inplace = True)
        cell_mean_list = []
        for loadno in load_sample.index:
            shift_blist = []
            for chrno in self.CHROMOSOME:
                #for chrno in ['chr20']:
                #chrno = 'chr20'
                read_file_path = Tads.load_path(self.config['sc_tad_path'] + '/', load_sample.loc[loadno, 'name'], res, chrno)
                try:
                    sctads = pd.read_csv(read_file_path, header=None, sep="\t", skiprows=1)
                except Exception:
                    print('Read ERROR of file ' + read_file_path)
                    continue
                sctads.columns = ['chr', 'startb', 'endb', 'name', 'value']
                #pdb.set_trace()
                sctads['midB'] = (sctads.startb + sctads.endb) // 2
                chr_ccd = common_bin[common_bin['chr']==chrno].copy()
                chr_ccd.reset_index(drop = True, inplace = True)
                sctads['inCad'] = 0
                for i in sctads.index:
                    withinCcd = chr_ccd[(sctads.loc[i, 'midB'] > chr_ccd.bin1) & (sctads.loc[i, 'midB'] < chr_ccd.bin2)].shape[0]
                    print('Cluster:', 'C' + str(clusterno+1), load_sample.loc[loadno, 'name'], chrno, i, withinCcd)
                    if withinCcd:
                        sctads.loc[i, 'inCad'] = 1
                if cads:
                    scintads = sctads[sctads.inCad==1].copy()
                else:
                    scintads = sctads[sctads.inCad==0].copy()
                scintads.reset_index(drop = True, inplace = True)
                if load_sample.loc[loadno, 'Sample'] == 'MCF7':
                    read_file_path = Tads.load_path(self.config['batch1_tad_path'], 'mcf7', res, chrno)
                    batch_tads = pd.read_csv(read_file_path, header=None, sep="\t", skiprows=1)
                elif load_sample.loc[loadno, 'Sample'] == 'MCF7M1':
                    read_file_path = Tads.load_path(self.config['batch2_tad_path'], 'mcf7m1', res, chrno)
                    batch_tads = pd.read_csv(read_file_path, header=None, sep="\t", skiprows=1)
                elif load_sample.loc[loadno, 'Sample'] == 'MCF7TR':
                    read_file_path = Tads.load_path(self.config['batch3_tad_path'], 'mcf7tr', res, chrno)
                    batch_tads = pd.read_csv(read_file_path, header=None, sep="\t", skiprows=1)
                else:
                    print('Reading TADs ERROR.')
                    exit
                batch_tads.columns = ['chr', 'startb', 'endb', 'name', 'value']
                batch_tads['midB'] = (batch_tads.startb + batch_tads.endb) // 2
                for j in scintads.index:
                    shift_boundary = np.min(np.abs(scintads.loc[j, 'midB'] - batch_tads.midB))
                    shift_blist.append(shift_boundary)
            cell_mean_list.append(np.mean(shift_blist))
        return cell_mean_list
    
    def shift_list(self, sample_name, pdCluster, cluster_list, cluster_num, res_list = [50000, 100000, 200000, 500000]):
        '''
        Get the list of shifted TADS boundaries of CADs and NADs.
        '''
        shift_cads_list = []
        shift_nads_list = []
        for res in res_list:
            shiftCads = []
            for i in range(cluster_num):
                common_bin_list = self.common_binset(sample_name, cluster_list, i)
                shiftCads.append(self.cell_shift(sample_name, pdCluster, i, common_bin_list, res=res, cads=True))
            shift_cads_list.append(shiftcads)
            shift_nads = []
            for i in range(cluster_num):
                common_bin_list = self.common_binset(sample_name, cluster_list, i)
                shift_nads.append(self.cell_shift(sample_name, pdCluster, i, common_bin_list, res=res, cads=False))
            shift_nads_list.append(shift_nads)
        return shift_cads_list, shift_nads_list
