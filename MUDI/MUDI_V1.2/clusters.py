#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
###Clusters and sub-clusters
###V9.20210322 for single cells MCF7, MCF7-TamR, MCF7-T1M
###Sub-clusters with all data (from batch 2, 3, 5 to batch 2, 6, 5)

import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from cluster import cpu
#from cluster import cputsne
from mpl_toolkits.mplot3d import Axes3D
import random

os.chdir('/data/yufan/schic/cluster/analysis06')
os.getcwd()

res = 1000000
allchrno = ["chr" + str(i+1) for i in range(22)] + ['chrX']

samplename = pd.read_csv('/data/yufan/schic/schicnames.txt', header=0, sep="\t")
###Change the batch 3 to batch 6
samplename['Round'] = [6 if i==3 else i for i in samplename.Round]

###countcutoff must be > 1, otherwise cpu.hicluster_cpu will be error
countcutoff = 2

#MCF7
groupname = samplename[samplename.Sample=='MCF7'].copy()
groupname.reset_index(drop = True, inplace = True)
group1 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group1.append(groupname.Name[i])

#MCF7-T1M
groupname = samplename[samplename.Sample=='MCF7-T1M'].copy()
groupname.reset_index(drop = True, inplace = True)
group2 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group2.append(groupname.Name[i])

#MCF7-TamR
groupname = samplename[samplename.Sample=='MCF7-TamR'].copy()
groupname.reset_index(drop = True, inplace = True)
group3 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group3.append(groupname.Name[i])

###K562
k562list = ['102_K562-B', '110_K562-B', '123_K562-B', '141_K562-B',
	'156_K562-A', '215_K562-A', '217_K562-B', '220_K562-B', '221_K562-A', 
	'222_K562-B', '223_K562-B', '224_K562-B', '225_K562-A', '226_K562-B',
	'227_K562-A', '228_K562-B', '229_K562-A', '230_K562-A', '231_K562-A',
	'232_K562-A', '233_K562-A', '234_K562-A', '235_K562-A', '236_K562-A',
	'237_K562-A', '238_K562-A', '239_K562-A', '240_K562-A', '241_K562-A',
	'242_K562-A', '54_K562-B', '58_K562-B', '72_K562-B', '79_K562-B']
###	'K562_bulkA', 'K562_bulkB'
group4 = []
for i in k562list:
	runpath = "/data/yufan/epigenetics/venus/yufan/singlecells/k562/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group4.append(i)

###Bing Ren WTC11C6
WTC11C6list = ['cell0' + str(i) if i < 10 else 'cell' + str(i) for i in range(1,95)]

group5 = []
for i in WTC11C6list:
	runpath = "/data/yufan/schic/public/ren4dn/hicpro/WTC11C6/output/hic_results/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = 'WTC11C6_'+ i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group5.append(i)

group5 = ['WTC11C6_' + i for i in group5]


###Bing Ren WTC11C6
WTC11C28list = ['cell0' + str(i) if i < 10 else 'cell' + str(i) for i in range(1,95)]

group6 = []
for i in WTC11C28list:
	runpath = "/data/yufan/schic/public/ren4dn/hicpro/WTC11C28/output/hic_results/matrix/"
	rundata = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + i + "/raw/" + str(res) + "/" + i + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = 'WTC11C28_'+ i + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group6.append(i)

group6 = ['WTC11C28_' + i for i in group6]

#network = group1 + group2 + group3
network = group1 + group2 + group3 + group4 + group5 + group6

label = network
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
chrom = [str(i+1) for i in range(22)] + ['X']
chromsize = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
nc = 2
ndim = 20

###clusterlabel: cluster labels identified by KMeans, rawpca: the matrix transformed by PCA
#clusterlabel, rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)
rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)

###make plot
pc1 = 0
pc2 = 1
pc3 = 2

cellnum1 = len(group1)
cellnum2 = len(group2)
cellnum3 = len(group3)
cellnum4 = len(group4)
cellnum5 = len(group5)
cellnum6 = len(group6)

#cellnum4 = len(mcf7tbatch2)
print(cellnum1)
print(cellnum2)
print(cellnum3)
print(cellnum4)
print(cellnum5)
print(cellnum6)

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1,pc1],rawpca[1][0:cellnum1,pc2],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2),pc1],rawpca[1][cellnum1:(cellnum1 + cellnum2),pc2],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc1],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc2],'r.')
p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc2],'y.')
p5, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc2],'c.')
p6, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc2],'k.')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper center')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
plt.xlabel('PC1', fontsize=18)
plt.ylabel('PC2', fontsize=18)
#plt.axis('square')
plt.show()

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1, pc2],rawpca[1][0:cellnum1, pc3],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2],rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3],'r.')
p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3],'y.')
p5, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5),pc3],'c.')
p6, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6),pc3],'k.')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper right')
plt.xlabel('PC2', fontsize=18)
plt.ylabel('PC3', fontsize=18)
#plt.axis('square')
plt.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
####separate the data point to three parts, marked them with different color
p1 = ax.scatter(rawpca[1][0:cellnum1, pc1], rawpca[1][0:cellnum1, pc2], rawpca[1][0:cellnum1, pc3], c='g')  ###draw the data points
p2 = ax.scatter(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc1], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3], c='b')  ###draw the data points
p3 = ax.scatter(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc1], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3], c='r')  ###draw the data points
p4 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3], c='y')  ###draw the data points
p5 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5), pc3], c='c')  ###draw the data points
p6 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5):(cellnum1 + cellnum2 + cellnum3 + cellnum4 + cellnum5 + cellnum6), pc3], c='k')  ###draw the data points

ax.set_zlabel('PC3')  ###axis
ax.set_ylabel('PC2')
ax.set_xlabel('PC1')
#ax.set_xticks(fontsize=10)
#ax.set_yticks(fontsize=10)
#ax.set_zticks(fontsize=10)
#Sometimes the default viewing angle is not optimal, 
#in which case we can use the view_init method to set the elevation and azimuthal angles. 
#In the following example, we'll use an elevation of 60 degrees (that is, 60 degrees above the x-y plane) 
#and an azimuth of 35 degrees (that is, rotated 35 degrees counter-clockwise about the z-axis):
#ax.view_init(60, 35)
ax.view_init(18, 161)
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
plt.legend([p1, p2, p3, p4, p5, p6], ['MCF7', 'MCF7M1', 'MCF7TR', 'K562', 'WTC11C6', 'WTC11C28'], loc='upper left', fontsize=8)

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
plt.show()

###############################################
###Make the sub-cluster

from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot
# 定义数据集
#X, _ = make_classification(n_samples=1000, n_features=2, n_informative=2, n_redundant=0, n_clusters_per_class=1, random_state=4)
# 定义模型
#X = rawpca[1][:,0:2]
X = rawpca[1][:,0:3]

#from sklearn.cluster import KMeans
#model = KMeans(n_clusters=7)

import warnings
warnings.filterwarnings("ignore")


#Calinski Harabaz Score
'''
from sklearn.metrics import calinski_harabaz_score
for i in range(12,1,-1):
	kmeans=GaussianMixture(n_components=i,random_state=123).fit(X)
	score=calinski_harabaz_score(X, kmeans.predict(X))
	print('The calinski_harabaz score for %d clustering is: %f'%(i,score))

randomseed = []
for i in range(1, 1001, 1):
	kmeans = GaussianMixture(n_components=7,random_state=i).fit(X)
	score = calinski_harabaz_score(X, kmeans.predict(X))
	randomseed.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

minseed = np.argmax(pd.Series(randomseed))
'''

#Silhouette Coefficient
from sklearn.metrics import silhouette_score
scorelist = []
for i in range(20,1,-1):
	randomscore = 0
	for j in range(1, 1001, 1):
		#print('i: ' + str(i) + ', j: ' + str(j))
		kmeans = GaussianMixture(n_components=i,random_state=j).fit(X)
		randomscore = randomscore + silhouette_score(X, kmeans.predict(X))
	scorelist.append(randomscore/1000)
	print('The calinski_harabaz score for %d clustering is: %f'%(i, randomscore/1000))

plt.plot(range(20, 1, -1), scorelist, color='black')
plt.xlim(20, 1)
plt.xticks(range(20, 1, -1), range(20, 1, -1), fontsize=13)
plt.xlabel('Cluster #', fontsize=15)
plt.ylabel('Silhouette Coefficient', fontsize=15)
plt.show()

scorelist = []
for i in range(1, 1001, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
	score = silhouette_score(X, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))

###4, 5, 6, 9, 10, 11, 13, 14 best 10
maxseed = 10
#高斯混合模型
model = GaussianMixture(n_components=9, random_state=maxseed)

#model = GaussianMixture(n_components=7, random_state=316)
#model = GaussianMixture(n_components=4, random_state=123)

################################
# 模型拟合
model.fit(X)
# 为每个示例分配一个集群
yhat = model.predict(X)
# 检索唯一群集
clusters = unique(yhat)
# 为每个群集的样本创建散点图
###PC1, PC2
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	pyplot.scatter(X[row_ix, 0], X[row_ix, 1])

# 绘制散点图
pyplot.xlabel('PC1', fontsize=18)
pyplot.ylabel('PC2', fontsize=18)
pyplot.show()

###PC2, PC3
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	pyplot.scatter(X[row_ix, 1], X[row_ix, 2])

# 绘制散点图
pyplot.xlabel('PC2', fontsize=18)
pyplot.ylabel('PC3', fontsize=18)
pyplot.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	ax.scatter(X[row_ix, 0], X[row_ix, 1], X[row_ix, 2])

# 绘制散点图
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.view_init(18, 161)
plt.show()

###Fowlkes-Mallows Score
labeltrue = [1] * cellnum1 + [2] * cellnum2 + [3] * cellnum3
from sklearn import metrics
scorelist = []
for i in range(20,1,-1):
	kmeans=GaussianMixture(n_components=i,random_state=316).fit(X)
#	score=metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The adjusted rand score for %d clustering is: %f'%(i,score))

scorelist = []
for i in range(1, 1000, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
#	score = metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))
###maxseed: 181

###
clusterlabel = ['MCF7'] * cellnum1 + ['MCF7-T1M'] * cellnum2 + ['MCF7-TamR'] * cellnum3
pdcluster = pd.DataFrame({'sample': network, 'label': clusterlabel, 'cluster': yhat, 'PC1': X[:, 0], 'PC2': X[:, 1], 'PC3': X[:, 2]})

pdcluster.groupby(['cluster', 'label']).size()
'''
cluster  label
0        MCF7-TamR     8
1        MCF7         39
2        MCF7-T1M     43
         MCF7-TamR    11
3        MCF7-T1M      1
         MCF7-TamR    14
4        MCF7         31
5        MCF7         15
         MCF7-T1M     10
         MCF7-TamR    16
6        MCF7-TamR    25
7        MCF7          2
         MCF7-TamR     8
8        MCF7-TamR     8
dtype: int64

'''


mcf7c5 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==5)]['sample']
mcf7t1mc5 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==5)]['sample']
mcf7tamrc5 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==5)]['sample']

mcf7t1mc2 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==2)]['sample']
mcf7tamrc2 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==2)]['sample']

mcf7c7 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==7)]['sample']
mcf7tamrc7 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==7)]['sample']

clusterc0 = pdcluster[pdcluster.cluster==0]['sample']
clusterc1 = pdcluster[pdcluster.cluster==1]['sample']
clusterc2 = pdcluster[pdcluster.cluster==2]['sample']
clusterc3 = pdcluster[pdcluster.cluster==3]['sample']
clusterc4 = pdcluster[pdcluster.cluster==4]['sample']
clusterc5 = pdcluster[pdcluster.cluster==5]['sample']
clusterc6 = pdcluster[pdcluster.cluster==6]['sample']
clusterc7 = pdcluster[pdcluster.cluster==7]['sample']
clusterc8 = pdcluster[pdcluster.cluster==8]['sample']

#common bins
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c5) + list(mcf7t1mc5) + list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7t1mc2) + list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc7))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c7) + list(mcf7tamrc7))].copy()
groupname = samplename[samplename.Name.isin(list(clusterc7) + list(clusterc8))].copy()
groupname.reset_index(drop = True, inplace = True)
setlist = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	csvmat['twobin'] = csvmat.bin1 * 10000+ csvmat.bin2
	setlist.append(set(csvmat.twobin))

commonset = setlist[0]			###get the commonset of all cells to commonset
for i in range(1,len(setlist)):
	commonset = commonset.intersection(setlist[i])

###calculate the common set
commonset = sorted(list(commonset))
###add the bins of chr22:2846-2897
#commonset.extend([i*10000+i for i in range(2846,2898)])
print(len(commonset))
###C5
###MCF7: 485
###MCF7-T1M: 173
###MCF7-TamR: 17
###All three: 8
###C2
###MCF7-T1M: 508
###MCF7-TamR: 1189
###All two: 268
###C7
###MCF7: 500
###MCF7-TamR: 6
###All two: 2




###cluster name and label number
###SC1: 3, SC2: 6, SC3: 1, SC4: 0, SC5: 4, SC6: 5, SC7: 2
pdcluster['sclabel'] = pdcluster.apply(lambda row : \
	'sc1' if row['cluster']==3 else \
	'sc2' if row['cluster']==6 else \
	'sc3' if row['cluster']==1 else \
	'sc4' if row['cluster']==0 else \
	'sc5' if row['cluster']==4 else \
	'sc6' if row['cluster']==5 else \
	'sc7' if row['cluster']==2 else \
	'ERROR', axis=1)

pdcluster.groupby(['sclabel', 'label']).size()
'''
sclabel  label
sc1      MCF7         27
sc2      MCF7         44
sc3      MCF7         15
         MCF7-T1M     10
         MCF7-TamR    16
sc4      MCF7-T1M     43
         MCF7-TamR     8
sc5      MCF7-TamR    28
sc6      MCF7-TamR    17
sc7      MCF7          1
         MCF7-T1M      1
         MCF7-TamR    21
dtype: int64
'''

pdcluster[pdcluster.sclabel=='sc3']










#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
###Clusters and sub-clusters
###V8.20210315 for single cells MCF7, MCF7-TamR, MCF7-T1M
###Sub-clusters with all data (from batch 2, 3, 5 to batch 2, 6, 5)
###Good version for clusters # 9

import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from cluster import cpu
#from cluster import cputsne
from mpl_toolkits.mplot3d import Axes3D
import random

os.chdir('/data/yufan/schic/cluster/analysis06')
os.getcwd()

res = 1000000
allchrno = ["chr" + str(i+1) for i in range(22)] + ['chrX']

samplename = pd.read_csv('/data/yufan/schic/schicnames.txt', header=0, sep="\t")
###Change the batch 3 to batch 6
samplename['Round'] = [6 if i==3 else i for i in samplename.Round]

###countcutoff must be > 1, otherwise cpu.hicluster_cpu will be error
countcutoff = 6

#MCF7
groupname = samplename[samplename.Sample=='MCF7'].copy()
groupname.reset_index(drop = True, inplace = True)
group1 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group1.append(groupname.Name[i])

#MCF7-T1M
groupname = samplename[samplename.Sample=='MCF7-T1M'].copy()
groupname.reset_index(drop = True, inplace = True)
group2 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group2.append(groupname.Name[i])

#MCF7-TamR
groupname = samplename[samplename.Sample=='MCF7-TamR'].copy()
groupname.reset_index(drop = True, inplace = True)
group3 = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	bedfile = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res) + "_abs.bed"
	print("Reading " + bedfile + "......")
	csvbed = pd.read_csv(bedfile, header=None, sep="\t")
	csvbed.columns = ["chr", "startb", "endb", "binnum"]
	chrcounts = []
	for chrno in allchrno:
		maxnum = max(csvbed[csvbed['chr']==chrno]['binnum'])
		minnum = min(csvbed[csvbed['chr']==chrno]['binnum'])
		cellchr = csvmat[(csvmat['bin1'] >= minnum) & (csvmat['bin1'] <= maxnum) & (csvmat['bin2'] >= minnum) & (csvmat['bin2'] <= maxnum)]
		cellchrnew = pd.concat([cellchr.iloc[:,0:2] - minnum, cellchr.iloc[:,2]], axis = 1)
		savefilename = groupname.Name[i] + "_" + chrno + ".txt"
		chrcounts.append(cellchrnew.shape[0])
		cellchrnew.to_csv(savefilename,sep='\t',header=False,index=False)
	addsample = True
	for j in chrcounts:
		if j < countcutoff:
			addsample = False
	if addsample:
		group3.append(groupname.Name[i])


network = group1 + group2 + group3
#network = group1 + group2 + group3 + group4

label = network
hg19dim = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560]
chrom = [str(i+1) for i in range(22)] + ['X']
chromsize = {chrom[i]:hg19dim[i] for i in range(len(chrom))}
nc = 2
ndim = 20

###clusterlabel: cluster labels identified by KMeans, rawpca: the matrix transformed by PCA
#clusterlabel, rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)
rawpca = cpu.hicluster_cpu(network, chromsize, nc=nc, res=res, ncpus=5)

###make plot
pc1 = 0
pc2 = 1
pc3 = 2

cellnum1 = len(group1)
cellnum2 = len(group2)
cellnum3 = len(group3)
#cellnum4 = len(group4)

#cellnum4 = len(mcf7tbatch2)
print(cellnum1)
print(cellnum2)
print(cellnum3)
#print(cellnum4)

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1,pc1],rawpca[1][0:cellnum1,pc2],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2),pc1],rawpca[1][cellnum1:(cellnum1 + cellnum2),pc2],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc1],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3),pc2],'r.')
#p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc1],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4),pc2],'y.')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower left')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
#plt.legend([p1, p2, p3, p4], ['MCF7', 'MCF7-T1M', 'MCF7-TamR', 'K562'], loc='upper left')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
plt.xlabel('PC1', fontsize=18)
plt.ylabel('PC2', fontsize=18)
#plt.axis('square')
plt.show()

plt.figure()
p1, = plt.plot(rawpca[1][0:cellnum1, pc2],rawpca[1][0:cellnum1, pc3],'g.')
p2, = plt.plot(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2],rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3],'b.')
p3, = plt.plot(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2],rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3],'r.')
#p4, = plt.plot(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2],rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3],'y.')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper left')
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='upper right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='lower right')
#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper center')
#plt.legend([p1, p2], ['MCF7', 'MCF7-TamR'], loc='upper right')
#plt.legend([p1, p2, p3, p4], ['MCF7', 'MCF7-T1M', 'MCF7-TamR', 'K562'], loc='upper right')
plt.xlabel('PC2', fontsize=18)
plt.ylabel('PC3', fontsize=18)
#plt.axis('square')
plt.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
####separate the data point to three parts, marked them with different color
p1 = ax.scatter(rawpca[1][0:cellnum1, pc1], rawpca[1][0:cellnum1, pc2], rawpca[1][0:cellnum1, pc3], c='g')  ###draw the data points
p2 = ax.scatter(rawpca[1][cellnum1:(cellnum1 + cellnum2), pc1], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc2], rawpca[1][cellnum1:(cellnum1 + cellnum2), pc3], c='b')  ###draw the data points
p3 = ax.scatter(rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc1], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc2], rawpca[1][(cellnum1 + cellnum2):(cellnum1 + cellnum2 + cellnum3), pc3], c='r')  ###draw the data points
#p4 = ax.scatter(rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc1], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc2], rawpca[1][(cellnum1 + cellnum2 + cellnum3):(cellnum1 + cellnum2 + cellnum3 + cellnum4), pc3], c='y')  ###draw the data points

ax.set_zlabel('PC3')  ###axis
ax.set_ylabel('PC2')
ax.set_xlabel('PC1')
#ax.set_xticks(fontsize=10)
#ax.set_yticks(fontsize=10)
#ax.set_zticks(fontsize=10)
#Sometimes the default viewing angle is not optimal, 
#in which case we can use the view_init method to set the elevation and azimuthal angles. 
#In the following example, we'll use an elevation of 60 degrees (that is, 60 degrees above the x-y plane) 
#and an azimuth of 35 degrees (that is, rotated 35 degrees counter-clockwise about the z-axis):
#ax.view_init(60, 35)
ax.view_init(18, 161)
plt.legend([p1, p2, p3], ['MCF7', 'MCF7M1', 'MCF7TR'], loc='upper left')
#plt.legend([p1, p2, p3, p4], ['MCF7', 'MCF7-T1M', 'MCF7-TamR', 'K562'], loc='upper left')

#plt.legend([p1, p2, p3], ['MCF7', 'MCF7-T1M', 'MCF7-TamR'], loc='upper right')
plt.show()

###############################################
###Make the sub-cluster

from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot
# 定义数据集
#X, _ = make_classification(n_samples=1000, n_features=2, n_informative=2, n_redundant=0, n_clusters_per_class=1, random_state=4)
# 定义模型
#X = rawpca[1][:,0:2]
X = rawpca[1][:,0:3]

#from sklearn.cluster import KMeans
#model = KMeans(n_clusters=7)

import warnings
warnings.filterwarnings("ignore")


#Calinski Harabaz Score
'''
from sklearn.metrics import calinski_harabaz_score
for i in range(12,1,-1):
	kmeans=GaussianMixture(n_components=i,random_state=123).fit(X)
	score=calinski_harabaz_score(X, kmeans.predict(X))
	print('The calinski_harabaz score for %d clustering is: %f'%(i,score))

randomseed = []
for i in range(1, 1001, 1):
	kmeans = GaussianMixture(n_components=7,random_state=i).fit(X)
	score = calinski_harabaz_score(X, kmeans.predict(X))
	randomseed.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

minseed = np.argmax(pd.Series(randomseed))
'''

#Silhouette Coefficient
from sklearn.metrics import silhouette_score
scorelist = []
for i in range(20,1,-1):
	randomscore = 0
	for j in range(1, 1001, 1):
		#print('i: ' + str(i) + ', j: ' + str(j))
		kmeans = GaussianMixture(n_components=i,random_state=j).fit(X)
		randomscore = randomscore + silhouette_score(X, kmeans.predict(X))
	scorelist.append(randomscore/1000)
	print('The calinski_harabaz score for %d clustering is: %f'%(i, randomscore/1000))

plt.plot(range(20, 1, -1), scorelist, color='black')
plt.xlim(20, 1)
plt.xticks(range(20, 1, -1), range(20, 1, -1), fontsize=13)
plt.xlabel('Cluster #', fontsize=15)
plt.ylabel('Silhouette Coefficient', fontsize=15)
plt.show()

scorelist = []
for i in range(1, 1001, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
	score = silhouette_score(X, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))

####################
###Fowlkes-Mallows Score
labeltrue = [1] * cellnum1 + [2] * cellnum2 + [3] * cellnum3
from sklearn import metrics
scorelist = []
for i in range(20,1,-1):
	kmeans=GaussianMixture(n_components=i,random_state=316).fit(X)
#	score=metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The adjusted rand score for %d clustering is: %f'%(i,score))

scorelist = []
for i in range(1, 500, 1):
	kmeans = GaussianMixture(n_components=9,random_state=i).fit(X)
#	score = metrics.adjusted_rand_score(labeltrue, kmeans.predict(X))
	score = metrics.fowlkes_mallows_score(labeltrue, kmeans.predict(X))
	scorelist.append(score)
	print('The calinski_harabaz score for random seed of %d is: %f'%(i,score))

maxseed = np.argmax(pd.Series(scorelist))
###maxseed: 181

###4, 5, 6, 9, 10, 11, 13, 14 best 10
maxseed = 10
#高斯混合模型
model = GaussianMixture(n_components=9, random_state=maxseed)

###Optimized by Fowlkes-Mallows Score
#Good: 12, 14, 17(!), 23(!), 24(!!), 29(!!!)
model = GaussianMixture(n_components=9, random_state=29)
#model = GaussianMixture(n_components=4, random_state=123)

################################
# 模型拟合
model.fit(X)
# 为每个示例分配一个集群
yhat = model.predict(X)
# 检索唯一群集
clusters = unique(yhat)
# 为每个群集的样本创建散点图
###PC1, PC2
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	pyplot.scatter(X[row_ix, 0], X[row_ix, 1])

# 绘制散点图
pyplot.xlabel('PC1', fontsize=18)
pyplot.ylabel('PC2', fontsize=18)
pyplot.show()

###PC2, PC3
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	pyplot.scatter(X[row_ix, 1], X[row_ix, 2])

# 绘制散点图
pyplot.xlabel('PC2', fontsize=18)
pyplot.ylabel('PC3', fontsize=18)
pyplot.show()

###3D figure
plt.figure()
ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
for cluster in clusters:
	# 获取此群集的示例的行索引
	row_ix = where(yhat == cluster)
	# 创建这些样本的散布
	ax.scatter(X[row_ix, 0], X[row_ix, 1], X[row_ix, 2])

# 绘制散点图
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.view_init(18, 161)
plt.show()

###
clusterlabel = ['MCF7'] * cellnum1 + ['MCF7-T1M'] * cellnum2 + ['MCF7-TamR'] * cellnum3
pdcluster = pd.DataFrame({'sample': network, 'label': clusterlabel, 'cluster': yhat, 'PC1': X[:, 0], 'PC2': X[:, 1], 'PC3': X[:, 2]})
pdcluster.to_csv('/data/yufan/schic/cluster/clusterlist001.txt',sep='\t',header=True,index=False)

pdcluster.groupby(['cluster', 'label']).size()
'''
cluster  label
0        MCF7-TamR     8
1        MCF7         39
2        MCF7-T1M     43
         MCF7-TamR    11
3        MCF7-T1M      1
         MCF7-TamR    14
4        MCF7         31
5        MCF7         15
         MCF7-T1M     10
         MCF7-TamR    16
6        MCF7-TamR    25
7        MCF7          2
         MCF7-TamR     8
8        MCF7-TamR     8
dtype: int64

'''


mcf7c5 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==5)]['sample']
mcf7t1mc5 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==5)]['sample']
mcf7tamrc5 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==5)]['sample']

mcf7t1mc2 = pdcluster[(pdcluster.label=='MCF7-T1M') & (pdcluster.cluster==2)]['sample']
mcf7tamrc2 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==2)]['sample']

mcf7c7 = pdcluster[(pdcluster.label=='MCF7') & (pdcluster.cluster==7)]['sample']
mcf7tamrc7 = pdcluster[(pdcluster.label=='MCF7-TamR') & (pdcluster.cluster==7)]['sample']

clusterc0 = pdcluster[pdcluster.cluster==0]['sample']
clusterc1 = pdcluster[pdcluster.cluster==1]['sample']
clusterc2 = pdcluster[pdcluster.cluster==2]['sample']
clusterc3 = pdcluster[pdcluster.cluster==3]['sample']
clusterc4 = pdcluster[pdcluster.cluster==4]['sample']
clusterc5 = pdcluster[pdcluster.cluster==5]['sample']
clusterc6 = pdcluster[pdcluster.cluster==6]['sample']
clusterc7 = pdcluster[pdcluster.cluster==7]['sample']
clusterc8 = pdcluster[pdcluster.cluster==8]['sample']

#common bins
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c5) + list(mcf7t1mc5) + list(mcf7tamrc5))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7t1mc2) + list(mcf7tamrc2))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7tamrc7))].copy()
#groupname = samplename[samplename.Name.isin(list(mcf7c7) + list(mcf7tamrc7))].copy()
#groupname = samplename[samplename.Name.isin(list(clusterc7) + list(clusterc8))].copy()
groupname = samplename[samplename.Name.isin(list(clusterc4))].copy()
#groupname = samplename[samplename.Name.isin(list(clusterc5))].copy()
groupname.reset_index(drop = True, inplace = True)
setlist = []
for i in groupname.index:
	runpath = "/data/yufan/schic/batch" + str(groupname.Round[i]) + "/hicpro2/output/hic_results/matrix/"
	rundata = runpath + groupname.Name[i] + "/raw/" + str(res) + "/" + groupname.Name[i] + "_" + str(res)+ ".matrix"
	print("Reading " + rundata + "......")
	csvmat = pd.read_csv(rundata, header=None, sep="\t")
	csvmat.columns = ["bin1", "bin2", "binvalue"]
	csvmat['twobin'] = csvmat.bin1 * 10000+ csvmat.bin2
	setlist.append(set(csvmat.twobin))

commonset = setlist[0]			###get the commonset of all cells to commonset
for i in range(1,len(setlist)):
	commonset = commonset.intersection(setlist[i])

###calculate the common set
commonset = sorted(list(commonset))
###add the bins of chr22:2846-2897
#commonset.extend([i*10000+i for i in range(2846,2898)])
print(len(commonset))
###C5
###MCF7: 485
###MCF7-T1M: 173
###MCF7-TamR: 17
###All three: 8
###C2
###MCF7-T1M: 508
###MCF7-TamR: 1189
###All two: 268
###C7
###MCF7: 500
###MCF7-TamR: 6
###All two: 2

####Gene on the loci
refseq = pd.read_csv('/data/yufan/biodb/refseq/refGene.txt', header=None, sep='\t')
refseq.columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
				'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score',
				'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']

refseqbed = refseq.loc[:,['name2', 'chrom', 'strand', 'txStart', 'txEnd']].copy()
loci1 = pd.Series(commonset)//10000
loci2 = pd.Series(commonset)%10000
for i in loci1.index:
	print(i)
	if loci1[i] == loci2[i]:
		print('True')
	else:
		print('False')

binbed = pd.read_csv('/data/yufan/schic/batch6/hicpro2/output/hic_results/matrix/TA01/raw/1000000/TA01_1000000_abs.bed', header=None, sep='\t')
binbed.columns = ['chr', 'bin1', 'bin2', 'pos']
commonbin = binbed[binbed.pos.isin(loci1)].copy()

import zhou
chr1 = list(refseqbed.chrom)
start1 = list(refseqbed.txStart)
end1 = list(refseqbed.txEnd)
chr2 = list(commonbin.chr)
start2 = list(commonbin.bin1)
end2 = list(commonbin.bin2)
target = zhou.overlap(chr1, start1, end1, chr2, start2, end2)
refseqbed['common'] = 0
refseqbed.loc[target[1],'common'] = 1

commongene = list(unique(list(refseqbed[refseqbed.common==1].name2)))
len(commongene)
#1907

cdmole = pd.read_csv('/data/yufan/biodb/hgnc/cd_group-471.csv', header=0, sep=',')
cdlist = list(cdmole['Approved symbol'])
set(commongene).intersection(set(cdlist))
sorted(list(set(commongene).intersection(set(cdlist))))
#['ABCB1', 'CD14', 'CD46', 'CD55', 'CR1', 'CR2', 'FCAMR', 'GP5', 'IL5RA', 'ITGB4', 'NCR2', 'PROM1', 'PTPRC', 'SDC2', 'TREM1']

surface = pd.read_csv('/data/yufan/biodb/surfaceproteins/surface.csv', header=None, sep='\t')
surface.columns = ['gene']
sorted(list(set(commongene).intersection(set(list(surface.gene)))))
#['ABCB1', 'ABCB4', 'ADAM22', 'ADAM23', 'ADORA3', 'AMIGO1', 'ANO1', 'ANO6', 'ATP13A3', 'BACE2', 'CA4', 'CD14', 'CD46', 'CD55', 'CDH10', 'CDH12', 'CDH13', 'CDH26', 'CDH4', 'CDH6', 'CDH9', 'CDHR3', 'CHL1', 'CHRM2', 'CLDN1', 'CLDN16', 'CLEC14A', 'CLEC5A', 'CNTN1', 'CNTN4', 'CNTN6', 'CNTNAP4', 'CPM', 'CR1', 'CR2', 'CRB2', 'CSF1', 'CSPG4', 'DCHS2', 'DCSTAMP', 'DDR2', 'DPP6', 'DSCAM', 'FCAMR', 'FSHR', 'FZD6', 'GABRB2', 'GJC3', 'GLDN', 'GLP1R', 'GP5', 'GPC2', 'GPR1', 'GPR22', 'GPR37', 'GPR61', 'GPR65', 'GRIA1', 'GRID2', 'GRIK2', 'GRM7', 'GRM8', 'HEG1', 'HTR1B', 'HTR1F', 'HTR2A', 'IL17RC', 'IL17RE', 'IL1RAP', 'IL5RA', 'ITGA9', 'ITGB4', 'ITGB5', 'KCNK5', 'KCNMB4', 'LAMP5', 'LHFPL4', 'LINGO4', 'LRFN5', 'LRP12', 'LRP1B', 'LRP8', 'LRRC15', 'LRRC4', 'LRRC4C', 'LRRN1', 'LRRN2', 'LRRN3', 'MC3R', 'MDGA2', 'MFAP3', 'MFSD6', 'MR1', 'MUC13', 'NCR2', 'NFASC', 'NIPAL2', 'NKAIN2', 'NKAIN3', 'NLGN1', 'NPSR1', 'NRCAM', 'NRXN3', 'NTNG1', 'OPRM1', 'OR10G2', 'OR10G3', 'OR2AE1', 'OR4D1', 'OR4D2', 'OR4E1', 'OR4E2', 'OR9A4', 'OXTR', 'PCDHA1', 'PCDHA10', 'PCDHA11', 'PCDHA12', 'PCDHA13', 'PCDHA2', 'PCDHA3', 'PCDHA4', 'PCDHA5', 'PCDHA6', 'PCDHA7', 'PCDHA8', 'PCDHA9', 'PCDHAC1', 'PCDHAC2', 'PCDHB1', 'PCDHB10', 'PCDHB11', 'PCDHB12', 'PCDHB13', 'PCDHB14', 'PCDHB15', 'PCDHB16', 'PCDHB2', 'PCDHB3', 'PCDHB4', 'PCDHB5', 'PCDHB6', 'PCDHB7', 'PCDHB8', 'PCDHB9', 'PCDHGA1', 'PCDHGA10', 'PCDHGA11', 'PCDHGA12', 'PCDHGA2', 'PCDHGA3', 'PCDHGA4', 'PCDHGA5', 'PCDHGA6', 'PCDHGA7', 'PCDHGA8', 'PCDHGA9', 'PCDHGB1', 'PCDHGB2', 'PCDHGB3', 'PCDHGB4', 'PCDHGB5', 'PCDHGB6', 'PCDHGB7', 'PCDHGC3', 'PCDHGC4', 'PCDHGC5', 'PIGR', 'PILRA', 'PILRB', 'PKHD1L1', 'PMEPA1', 'PROM1', 'PRRT3', 'PTGER4', 'PTPRB', 'PTPRC', 'PTPRG', 'PTPRM', 'PTPRZ1', 'RNF43', 'RNFT1', 'ROBO2', 'SCN3A', 'SDC2', 'SDK2', 'SEMA6C', 'SEMA6D', 'SLC12A1', 'SLC12A8', 'SLC16A4', 'SLC16A5', 'SLC1A7', 'SLC24A5', 'SLC26A3', 'SLC26A4', 'SLC26A5', 'SLC28A1', 'SLC2A13', 'SLC2A9', 'SLC32A1', 'SLC38A11', 'SLC5A7', 'SLC5A8', 'SLC6A17', 'SLC6A5', 'SLCO1B3', 'SLCO1C1', 'SPPL2A', 'SSTR1', 'SSTR2', 'STEAP4', 'SUSD3', 'SYNPR', 'SYPL1', 'TAS2R3', 'TAS2R38', 'TAS2R4', 'TMEM178B', 'TMEM30A', 'TMEM63C', 'TMEM67', 'TMX4', 'TREM1', 'TREM2', 'TREML2', 'TRHR', 'TSHR', 'TSPAN2', 'TSPAN5', 'USH2A']
len(sorted(list(set(commongene).intersection(set(list(surface.gene))))))
#231

rnaseq = pd.read_csv('/data/yufan/epigenetics/backup/yufan/data20161003_mcf7_t47d_tamoxifen_resistant_rnaseq/results_tamr_mcf7_rnaseq/diff_out/gene_exp.diff', header=0, sep='\t')
tamrdeg = rnaseq[(rnaseq['p_value'] < 0.01) & (rnaseq['log2(fold_change)'] > 0.6)].copy()

targetdeg = tamrdeg[tamrdeg.gene.isin(commongene)].copy()
targetdeg.reset_index(drop = True, inplace = True)

sorted(list(targetdeg.gene))
###Cluster 4
#['AASS', 'ACOX1', 'ADAMTS7', 'ADAMTS9', 'ADGRA3', 'ADH5', 'AKAP6', 'AMPD1', 'ARL14EP', 'ATP8B4', 'BAALC-AS1', 'BACE2', 'BACH2', 'BAG2', 'BCAS1', 'BCAS2', 'BMP7', 'BRWD1', 'C1orf116', 'C3orf14', 'C3orf35', 'CADPS2', 'CCSER1', 'CCT2', 'CGN', 'CMC1', 'CNKSR3', 'CNOT2', 'COG3', 'COX7A2', 'CRELD1', 'CTHRC1', 'CTSH', 'CUEDC1', 'DDX27', 'DECR1', 'DENND2C', 'DPM1', 'DPY19L1', 'EAPP', 'EEF1B2', 'EGLN3', 'EIF3E', 'ERG', 'ETS2', 'FAM200A', 'FAM219B', 'FAM43A', 'FEZF1', 'FEZF1-AS1', 'FOXA1', 'FOXP4', 'GAL3ST4', 'GALNT10', 'GDPD1', 'GLS', 'GOLGA6L3', 'GPR37', 'HEY2', 'HIF1A', 'HIF1A-AS2', 'HOXD8', 'INPP1', 'JADE1', 'KCNC4', 'KCND3', 'KCNG1', 'KCNMB4', 'LAPTM4B', 'LDB2', 'LINC00643', 'LINC00994', 'LOC101929709', 'LRP12', 'LRRC6', 'MBIP', 'MDGA2', 'MECOM', 'MFSD6', 'MIPOL1', 'MIR4712', 'MKRN2OS', 'MPHOSPH6', 'MR1', 'MRPL9', 'MYRFL', 'NAB1', 'NAV2', 'NAV2-AS2', 'NDRG1', 'NIPAL2', 'NPAS3', 'NRCAM', 'OXR1', 'OXTR', 'PABPC1', 'PCDHB10', 'PCDHB11', 'PCDHB12', 'PCDHB13', 'PCDHB2', 'PCDHB3', 'PDE4D', 'PDE8A', 'PFDN4', 'PFKFB2', 'PGRMC2', 'PLCB4', 'PLRG1', 'PPARG', 'PRICKLE1', 'PRICKLE2', 'PRICKLE2-AS1', 'PRKAA1', 'PRKCH', 'PRR11', 'PSMB4', 'PSMD6', 'PTCSC3', 'PTGER4', 'PTPN21', 'PTPRB', 'PTPRG', 'RAB5A', 'RALGAPA2', 'RAP1GDS1', 'RFTN1', 'RIMS2', 'RNF19A', 'RNFT1', 'RPL30', 'RPS6KB1', 'RRBP1', 'S100A10', 'SALL2', 'SATB1', 'SCAMP5', 'SCP2', 'SDK2', 'SELENBP1', 'SEMA3E', 'SEMA6C', 'SLC12A8', 'SLC27A2', 'SLC2A13', 'SLC39A11', 'SNTB1', 'SNTN', 'SNX27', 'SPAG1', 'SPOCK3', 'SPTB', 'SPTSSA', 'SSTR2', 'ST3GAL1', 'STAU1', 'SULF2', 'SUPT4H1', 'SYNDIG1', 'SYNPR-AS1', 'SYT16', 'TAPT1', 'TATDN1', 'TFEB', 'THOC7', 'THUMPD3-AS1', 'TMEM178B', 'TMEM30A', 'TMEM40', 'TMEM63C', 'TMTC2', 'TMX4', 'TOM1L1', 'TRAM1L1', 'TRIB1', 'TRIM33', 'TSPAN5', 'TTC33', 'TTC6', 'UQCRB', 'VASH1', 'VEZF1', 'VMP1', 'WDR1', 'WEE2-AS1', 'WNK2', 'YPEL2', 'ZCWPW2', 'ZFAS1', 'ZNFX1', 'ZSCAN2', 'ZSWIM6']

sorted(list(set(list(targetdeg.gene)).intersection(set(cdlist))))
#[]

sorted(list(set(list(targetdeg.gene)).intersection(set(list(surface.gene)))))
#['BACE2', 'GPR37', 'KCNMB4', 'LRP12', 'MDGA2', 'MFSD6', 'MR1', 'NIPAL2', 'NRCAM', 'OXTR', 'PCDHB10', 'PCDHB11', 'PCDHB12', 'PCDHB13', 'PCDHB2', 'PCDHB3', 'PTGER4', 'PTPRB', 'PTPRG', 'RNFT1', 'SDK2', 'SEMA6C', 'SLC12A8', 'SLC2A13', 'SSTR2', 'TMEM178B', 'TMEM30A', 'TMEM63C', 'TMX4', 'TSPAN5']
len(sorted(list(set(list(targetdeg.gene)).intersection(set(list(surface.gene))))))
#30

###Cluster 5
#['ADAMTS9', 'AMPD1', 'BCAS1', 'BCAS2', 'C3orf14', 'DENND2C', 'LINC00994', 'PFDN4', 'PRICKLE2', 'PRICKLE2-AS1', 'PSMD6', 'PTPRG', 'RNFT1', 'RPS6KB1', 'SULF2', 'TRIM33']

for i in list(targetdeg.gene):
	print(i)

'''
targetdeg.sort_values(by='gene', ascending=True, inplace=True)
targetdeg.gene.to_csv('tamrpdldeg.csv', index=False, header=False, sep='\t')
'''

rnaseq['gseascore'] = rnaseq['log2(fold_change)'] / rnaseq['p_value']
rnaseq.sort_values(by='gseascore', ascending=False, inplace=True)
commonexp = rnaseq[rnaseq.gene.isin(commongene)].copy()
savedeg = commonexp.loc[:, ['gene', 'gseascore']].copy()
savedeg.replace([np.inf, -np.inf], np.nan).dropna(subset=['gseascore'], how='all').to_csv('/data/yufan/schic/cluster/gsea_common_gene_cluster4.rnk', index=False, header=False, sep='\t')


















###cluster name and label number
###SC1: 3, SC2: 6, SC3: 1, SC4: 0, SC5: 4, SC6: 5, SC7: 2
pdcluster['sclabel'] = pdcluster.apply(lambda row : \
	'sc1' if row['cluster']==3 else \
	'sc2' if row['cluster']==6 else \
	'sc3' if row['cluster']==1 else \
	'sc4' if row['cluster']==0 else \
	'sc5' if row['cluster']==4 else \
	'sc6' if row['cluster']==5 else \
	'sc7' if row['cluster']==2 else \
	'ERROR', axis=1)

pdcluster.groupby(['sclabel', 'label']).size()
'''
sclabel  label
sc1      MCF7         27
sc2      MCF7         44
sc3      MCF7         15
         MCF7-T1M     10
         MCF7-TamR    16
sc4      MCF7-T1M     43
         MCF7-TamR     8
sc5      MCF7-TamR    28
sc6      MCF7-TamR    17
sc7      MCF7          1
         MCF7-T1M      1
         MCF7-TamR    21
dtype: int64
'''

pdcluster[pdcluster.sclabel=='sc3']

