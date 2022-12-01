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
        
    def getBinBed(self, bedPath = './data/samples/cell_1000000_abs.bed'):
        '''
        self.binBed: list of genome loci for bins.
        '''
        self.binBed = pd.read_csv(bedPath, header=None, sep='\t')
        self.binBed.columns = ['chr', 'bin1', 'bin2', 'pos']
    
    def cellLabel(self, labelPath = './data/samples/cellLabel.txt'):
        '''
        Load the file of cell labeling.
        self.sampleName: list of sample
        '''
        self.sampleName = pd.read_csv(labelPath, header=0, sep="\t")
    
    def sampleCount(self, sampleList=[], labelList=[]):
        '''
        Count the samples and make the plot.
        '''
        if not sampleList:
            sampleList = list(set(self.sampleName['sample']))
        sampleCount = []
        sampleCount.append(self.sampleName[self.sampleName['sample'] == sampleList[0]].shape[0])
        sampleCount.append(self.sampleName[self.sampleName['sample'] == sampleList[1]].shape[0])
        sampleCount.append(self.sampleName[self.sampleName['sample'] == sampleList[2]].shape[0])
        ###Plot the sample count of single cell Hi-C
        if not labelList:
            labelList = sampleList
        rects1 = plt.bar(labelList, sampleCount, color = ['#88224488', '#22884488', '#22448888'])
        
        def autolabel(rects):
            for rect in rects:
                height = rect.get_height()
                #plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, '%s' % float(height))
                plt.text(rect.get_x()+rect.get_width() / 2 - 0.1, 1.03*height, str(height))
        
        autolabel(rects1)
        plt.ylim(0, 125)
        plt.ylabel('# of cells', fontsize=15)
        plt.show()
    
    def cellCluster(self, clusterPath = './data/samples/clusterList.txt'):
        '''
        Load cell cluster file.
        self.pdCluster: list of clusters.
        '''
        self.pdCluster = pd.read_csv(clusterPath, sep='\t', header=0)
        self.pdCluster['label'] = ['C' + str(i+1) for i in self.pdCluster.cluster]
        self.clusterList = []
        self.clusterNum = len(set(self.pdCluster.cluster))
        for i in range(self.clusterNum):
            self.clusterList.append(list(self.pdCluster[self.pdCluster.cluster==i]['sample']))
    
    def commonBinNum(self, sampleName, clusterList, binValueCutOff = 0):
        '''
        Return the number of Commonly Associating Domain(CADs) of the cluster.
        '''
        groupName = sampleName[sampleName.name.isin(clusterList)].copy()
        groupName.reset_index(drop = True, inplace = True)
        setList = []
        #res = 2000000
        for i in groupName.index:
            runData = groupName.loc[i, 'path']
            print("Reading " + runData + "......")
            csvMat = pd.read_csv(runData, header=None, sep="\t")
            csvMat.columns = ["bin1", "bin2", "binValue"]
            csvMat = csvMat[csvMat.binValue > binValueCutOff]
            csvMat['twoBin'] = csvMat.bin1 * 10000+ csvMat.bin2
            setList.append(set(csvMat.twoBin))
        commonSet = setList[0]			###get the commonset of all cells to commonset
        for i in range(1,len(setList)):
            commonSet = commonSet.intersection(setList[i])
        ###calculate the common set
        commonSet = sorted(list(commonSet))
        ###add the bins of chr22:2846-2897
        #commonset.extend([i*10000+i for i in range(2846,2898)])
        #print(len(commonset))
        ###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
        ###
        return len(commonSet)
    
    def commonBinPerCell(self):
        '''
        CADs number per cell.
        '''
        perCell = []
        for i in range(self.clusterNum):
            cellNumber = self.pdCluster[self.pdCluster.cluster==i].shape[0]
            perCell.append(self.commonBinNum(self.sampleName, self.clusterList[i], 1000000) / cellNumber)
        return perCell
    
    def commonBinSet(self, sampleName, clusterList, binValueCutOff = 0):
        '''
        The set of CADs.
        '''
        groupName = sampleName[sampleName.name.isin(clusterList)].copy()
        groupName.reset_index(drop = True, inplace = True)
        setList = []
        #res = 2000000
        for i in groupName.index:
            runData = groupName.loc[i, 'path']
            print("Reading " + runData + "......")
            csvMat = pd.read_csv(runData, header=None, sep="\t")
            csvMat.columns = ["bin1", "bin2", "binValue"]
            csvMat = csvMat[csvMat.binValue > binValueCutOff]
            csvMat['twoBin'] = csvMat.bin1 * 10000+ csvMat.bin2
            setList.append(set(csvMat.twoBin))
        commonSet = setList[0]			###get the commonset of all cells to commonset
        for i in range(1,len(setList)):
            commonSet = commonSet.intersection(setList[i])
        ###calculate the common set
        commonSet = sorted(list(commonSet))
        ###add the bins of chr22:2846-2897
        #commonset.extend([i*10000+i for i in range(2846,2898)])
        #print(len(commonset))
        ###CL1: 295, CL2: 2, CL3: 669 , CL4: 1, CL5: 268 , CL6: 8, CL7: 3, CL8: 2, CL9: 321
        ###
        return commonSet
    
    def cadScore(self, sampleName, clusterList, binValueCutOff = 0):
        '''
        Make the score of CADs.
        self.binList: list of interaction frequencies of each cell in a cluster.
        self.binScore: data frame of interaction frequencies, common bins (# of cells) of a cluster.
        '''
        groupName = sampleName[sampleName.name.isin(clusterList)].copy()
        groupName.reset_index(drop = True, inplace = True)
        ###Prepare to get the list of clusters and store in the list of self.binlist
        self.binList = []
        for i in groupName.index:
            runData = groupName.loc[i, 'path']
            print("Reading " + runData + "......")
            csvMat = pd.read_csv(runData, header=None, sep="\t")
            csvMat.columns = ["bin1", "bin2", "binValue"]
            csvMat = csvMat[csvMat.binValue > binValueCutOff]
            csvMat['twoBin'] = csvMat.bin1 * pow(10, len(str(np.max(csvMat.bin2)))) + csvMat.bin2
            self.binList.append(csvMat)
        ###Generate the cluster bin with the number of cells in column of cell
        clusterBin = self.binList[0]
        clusterBin['cell'] = 1
        for i in range(len(self.binList))[1:]:
            print(i+1, '/', len(self.binList))
            nextBin = self.binList[i]
            nextBin['cell'] = 1
            commonList = list(set(clusterBin.twoBin).intersection(set(nextBin.twoBin)))
            firstDiff = clusterBin[~clusterBin.twoBin.isin(commonList)].copy()
            nextDiff = nextBin[~nextBin.twoBin.isin(commonList)].copy()
            ###reset the cell number of both common bin
            bothCommon = clusterBin[clusterBin.twoBin.isin(commonList)].copy()
            bothCommon['cell'] = bothCommon['cell'] + 1
            ###reset the bin value
            bothCommon2 = nextBin[nextBin.twoBin.isin(commonList)].copy()
            bothCommon.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            bothCommon2.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            bothCommon['binValue2'] = list(bothCommon2['binValue'])
            bothCommon['binValue'] = bothCommon.apply(lambda row : max(row['binValue'], row['binValue2']), axis=1)
            bothCommon.drop('binValue2', axis=1, inplace=True)
            ###combine three df
            clusterBin = pd.concat([firstDiff, nextDiff, bothCommon], axis = 0)
            clusterBin.sort_values(by=['bin1', 'bin2'], ascending=True, inplace=True)
            clusterBin.reset_index(drop = True, inplace = True)
        ###save the score to binscore
        self.binScore = clusterBin.copy()
    
    def cadsLoci(self, bedDf, scoreDf):
        '''
        Get the genes of CADs.
        beddf: self.binbed
        scoredf: self.binscore
        '''
        ###get bins within binbed
        binList = list(self.binBed.pos)
        scoreDf = scoreDf[scoreDf.bin1.isin(binList) & scoreDf.bin2.isin(binList)].copy()
        scoreDf.reset_index(drop = True, inplace = True)
        ###get the genome loci
        scoreDf['bin1Chr'] = ''
        scoreDf['bin1Start'] = 0
        scoreDf['bin1End'] = 0
        scoreDf['bin2Chr'] = ''
        scoreDf['bin2Start'] = 0
        scoreDf['bin2End'] = 0
        for i in scoreDf.index:
            print(i+1, '/', scoreDf.shape[0])
            scoreDf.loc[i, 'bin1Chr']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin1']]['chr'])[0]
            scoreDf.loc[i, 'bin1Start']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin1']]['bin1'])[0]
            scoreDf.loc[i, 'bin1End']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin1']]['bin2'])[0]
            scoreDf.loc[i, 'bin2Chr']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin2']]['chr'])[0]
            scoreDf.loc[i, 'bin2Start']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin2']]['bin1'])[0]
            scoreDf.loc[i, 'bin2End']= list(bedDf[bedDf.pos==scoreDf.loc[i, 'bin2']]['bin2'])[0]
        return scoreDf
    
    def clusterGene(self, expandedBin, mode='max'):
        '''
        Make genes frequencies of clusters from expanded bins which have genome loci.
        '''
        geneCell = self.refSeq.copy()
        ###initiate cell number and interaction frequency
        geneCell['cell'] = 0
        geneCell['binValue'] = 0
        for i in geneCell.index:
            ### [a.start,a.end] and [b.start,b.end] overlap?
            ### if(A.start <= B.end && B.start <= A.end) overlap
            ### if(A.start > B.end || B.start > A.end)   No overlap
            geneBin = expandedBin[((geneCell.loc[i, 'chrNo'] == expandedBin.bin1Chr) & (geneCell.loc[i, 'start'] <= expandedBin.bin1End) & (expandedBin.bin1Start <= geneCell.loc[i, 'end'])) |
                                  ((geneCell.loc[i, 'chrNo'] == expandedBin.bin2Chr) & (geneCell.loc[i, 'start'] <= expandedBin.bin2End) & (expandedBin.bin2Start <= geneCell.loc[i, 'end']))].copy()
            print(i+1, '/', geneCell.shape[0], geneCell.loc[i, 'gene'])
            #print(genebin)
            if geneBin.shape[0]!=0:
                if mode=='max':
                    geneCell.loc[i, 'cell'] = np.max(geneBin.cell)
                    geneCell.loc[i, 'binValue'] = np.max(geneBin.binValue)
                elif mode=='mean':
                    geneCell.loc[i, 'cell'] = np.mean(geneBin.cell)
                    geneCell.loc[i, 'binValue'] = np.mean(geneBin.binValue)
                elif mode=='min':
                    geneCell.loc[i, 'cell'] = np.min(geneBin.cell)
                    geneCell.loc[i, 'binValue'] = np.min(geneBin.binValue)
                else:
                    print("Wrong mode.")
                ###Normalize the binvalue
                geneCell.loc[i, 'binValue'] = np.log2(geneCell.loc[i, 'binValue'])
        geneCell.sort_values(by=['gene', 'cell'], ascending=False, inplace=True)
        geneCell = geneCell.drop_duplicates(subset='gene', keep='first', inplace=False).copy()
        geneCell.sort_values(by='gene', ascending=True, inplace=True)
        geneCell.reset_index(drop = True, inplace = True)
        return geneCell
    
    def clusterGeneList(self, cutoffList = []):
        '''
        Return the list of genes in CADs of each cluster to:
        self.geneList
        '''
        # reset cutoffList as [0, 0, 0, ...]
        if not cutoffList:
            cutoffList = [0] * len(self.clusterList)
        ### get list of genes in CADs
        self.geneList = []
        for i in range(len(self.clusterList)):
            self.cadScore(sampleName = self.sampleName, clusterList = self.clusterList[i], binValueCutOff = cutoffList[i])
            cadsBin = self.binScore[self.binScore.cell==len(self.clusterList[i])].copy()
            cadsBin.reset_index(drop = True, inplace = True)
            cadsBin = self.cadsLoci(self.binBed, cadsBin)
            geneCell = self.clusterGene(cadsBin, mode='max')
            self.geneList.append(geneCell[geneCell.cell==len(self.binList)].copy())
    
    def integrationScore(self):
        '''
        Return the gene list with integration score of each cluster
        allPdList[i][j]: i: cluster of scHi-C, j: cluster of scRNA-seq
        '''
        allPdList = []
        for i in range(len(self.clusterList)):
            clusterGeneCell = self.geneList[i]
            genePdList = []
            for j in range(np.max(self.deg.clusterDEG.cluster)+1):
                #pdb.set_trace()
                clusterDeg2 = self.deg.clusterDEG[self.deg.clusterDEG.cluster==j].copy()
                overlappedGene = set(clusterGeneCell.gene).intersection(set(clusterDeg2.gene))
                geneListSorted = sorted(list(overlappedGene))
                ###integration score
                listFg = clusterGeneCell[clusterGeneCell.gene.isin(overlappedGene)].sort_values(by='gene', ascending=True, inplace=False).binValue.copy()
                listEg = clusterDeg2[clusterDeg2.gene.isin(overlappedGene)].sort_values(by='gene', ascending=True, inplace=False).avg_log2FC.copy()
                listFg.reset_index(drop = True, inplace = True)
                listEg.reset_index(drop = True, inplace = True)
                integrationGene = listFg * listEg / self.deg.totalDEG[j] / self.deg.clusterNumber.numNorm[j] * np.sum(self.deg.totalDEG)
                if str(integrationGene)=='nan':
                    integrationGene = 0
                genePd = pd.DataFrame({'gene': geneListSorted, 'integration': integrationGene})
                genePd.sort_values(by='integration', ascending=False, inplace=True)
                genePdList.append(genePd)
            allPdList.append(genePdList)
        return allPdList
    
    def commonBinBed(self, sampleName, clusterList):
        '''
        Return the bed format of CADs.
        chr: chromosome, bin1: start, bin2: end, pos: bin number.
        self.binBed: bed file to parse genome loci of the bin.
        '''
        commonSet = self.commonBinSet(sampleName = sampleName, clusterList = clusterList)
        ###loci of CADs
        loci1 = pd.Series(commonSet)//10000
        loci2 = pd.Series(commonSet)%10000
        for i in loci1.index:
            print(i)
            if loci1[i] == loci2[i]:
                print('True')
            else:
                print('False')
        commonBinBed = self.binBed[self.binBed.pos.isin(loci1)].copy()
        return commonBinBed
    
    def cellShift(self, sampleName, pdCluster, clusterNo, commonBin, scTadPath = '/data/insulation', batch1TadPath = '/data/insulation/mcf7/', batch2TadPath = '/data/insulation/mcf7m1/', batch3TadPath = '/data/insulation/mcf7tr/', res=100000, cads=True):
        '''
        Get the shifted TADs within CADs.
        clusterno: Number of cluster
        commonBin: Common bin of cluster (CADs)
        res: Resolution for TADs
        cads: True, CADs, False: NADs 
        '''
        loadSample = sampleName[sampleName.name.isin(pdCluster[pdCluster.cluster==clusterNo]['sample'])].copy()
        loadSample.reset_index(drop = True, inplace = True)
        cellMeanList = []
        for loadNo in loadSample.index:
            shiftBList = []
            for chrNo in self.allChrNo:
                #for chrno in ['chr20']:
                #chrno = 'chr20'
                readFilePath = TAD.loadPath(scTadPath + '/', loadSample.loc[loadno, 'name'], res, chrNo)
                try:
                    scTADs = pd.read_csv(readFilePath, header=None, sep="\t", skiprows=1)
                except Exception:
                    print('Read ERROR of file ' + readFilePath)
                    continue
                scTADs.columns = ['chr', 'startB', 'endB', 'name', 'value']
                #pdb.set_trace()
                scTADs['midB'] = (scTADs.startB + scTADs.endB) // 2
                chrCcd = commonBin[commonBin['chr']==chrNo].copy()
                chrCcd.reset_index(drop = True, inplace = True)
                scTADs['inCad'] = 0
                for i in scTADs.index:
                    withinCcd = chrCcd[(scTADs.loc[i, 'midB'] > chrCcd.bin1) & (scTADs.loc[i, 'midB'] < chrCcd.bin2)].shape[0]
                    print('Cluster:', 'C' + str(clusterNo+1), loadSample.loc[loadNo, 'name'], chrNo, i, withinCcd)
                    if withinCcd:
                        scTADs.loc[i, 'inCad'] = 1
                if cads:
                    scInTads = scTADs[scTADs.inCad==1].copy()
                else:
                    scInTads = scTADs[scTADs.inCad==0].copy()
                scInTads.reset_index(drop = True, inplace = True)
                if loadSample.loc[loadNo, 'Sample'] == 'MCF7':
                    readFilePath = TAD.loadPath(batch1TadPath, 'mcf7', res, chrNo)
                    batchTADs = pd.read_csv(readFilePath, header=None, sep="\t", skiprows=1)
                elif loadSample.loc[loadNo, 'Sample'] == 'MCF7M1':
                    readFilePath = TAD.loadPath(batch2TadPath, 'mcf7m1', res, chrNo)
                    batchTADs = pd.read_csv(readFilePath, header=None, sep="\t", skiprows=1)
                elif loadSample.loc[loadNo, 'Sample'] == 'MCF7TR':
                    readFilePath = TAD.loadPath(batch3TadPath, 'mcf7tr', res, chrNo)
                    batchTADs = pd.read_csv(readFilePath, header=None, sep="\t", skiprows=1)
                else:
                    print('Reading TADs ERROR.')
                    exit
                batchTADs.columns = ['chr', 'startB', 'endB', 'name', 'value']
                batchTADs['midB'] = (batchTADs.startB + batchTADs.endB) // 2
                for j in scInTads.index:
                    shiftBoundary = np.min(np.abs(scInTads.loc[j, 'midB'] - batchTADs.midB))
                    shiftBList.append(shiftBoundary)
            cellMeanList.append(np.mean(shiftBList))
        return cellMeanList
    
    def shiftList(self, sampleName, pdCluster, clusterList, clusterNum, resList = [50000, 100000, 200000, 500000]):
        '''
        Get the list of shifted TADS boundaries of CADs and NADs.
        '''
        shiftCadsList = []
        shiftNadsList = []
        for res in resList:
            shiftCads = []
            for i in range(clusterNum):
                commonBinList = self.commonBin(sampleName, clusterList, i)
                shiftCads.append(self.cellshift(sampleName, pdCluster, i, commonBinList, res=res, cads=True))
            shiftCadsList.append(shiftcads)
            shiftNads = []
            for i in range(clusterNum):
                commonBinList = self.commonBin(sampleName, clusterList, i)
                shiftNads.append(self.cellShift(sampleName, pdCluster, i, commonBinList, res=res, cads=False))
            shiftNadsList.append(shiftNads)
        return shiftCadsList, shiftNadsList

def test(shiftRun=False):
    singleCell = scHiC()
    singleCell.cellLabel()
    singleCell.sampleCount(sampleList=['MCF7', 'MCF7M1', 'MCF7TR'], labelList=['MCF7', 'MCF7M1', 'MCF7TR'])
    singleCell.cellCluster()
    singleCell.commonBinNum(singleCell.sampleName, singleCell.clusterList[0])
    print([round(i, 2) for i in singleCell.commonBinPerCell()])
    singleCell.commonBinSet(singleCell.sampleName, singleCell.clusterList[0])
    singleCell.commonBin(singleCell.sampleName, singleCell.clusterList, 0)
    commonBinBed = singleCell.commonBin(singleCell.sampleName, singleCell.clusterList, 0)
    if shiftRun:
        singleCell.cellShift(singleCell.sampleName, singleCell.pdCluster, 0, commonBinBed, res=100000, cads=True)
        shiftCads, shiftNads = singleCell.shiftList(singleCell.sampleName, singleCell.pdCluster, singleCell.clusterList, singleCell.clusterNum, resList = [50000, 100000, 200000, 500000])

