### modified from scHiCluster
### please refer to: https://github.com/zhoujt1994/scHiCluster

import os
import sys
import time
import numpy as np
#import torch
#import torch.nn.functional as F
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from scipy.stats import chi2_contingency
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics.cluster import adjusted_rand_score as ARI
#import pandas as pd
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

def neighbor_ave_cpu(A, pad):
	if pad==0:
		return A
	ngene, _ = A.shape
	ll = pad * 2 + 1
	B, C, D, E = [np.zeros((ngene + ll, ngene + ll)) for i in range(4)]
	B[(pad + 1):(pad + ngene + 1), (pad + 1):(pad + ngene + 1)] = A[:]
	F = B.cumsum(axis = 0).cumsum(axis = 1)
	C[ll :, ll:] = F[:-ll, :-ll]
	D[ll:, :] = F[:-ll, :]
	E[:, ll:] = F[:, :-ll]
	return (np.around(F + C - D - E, decimals=8)[ll:, ll:] / float(ll * ll))

def random_walk_cpu(A, rp):
	ngene, _ = A.shape
	A = A - np.diag(np.diag(A))
	A = A + np.diag(np.sum(A, axis=0) == 0)
	P = np.divide(A, np.sum(A, axis = 0))
	Q = np.eye(ngene)
	I = np.eye(ngene)
	for i in range(30):
		Q_new = (1 - rp) * I + rp * np.dot(Q, P)
		delta = np.linalg.norm(Q - Q_new)
		Q = Q_new.copy()
		if delta < 1e-6:
			break
	return Q

def impute_cpu(args):
	cell, c, ngene, pad, rp = args
	D = np.loadtxt(cell + '_chr' + c + '.txt')
	A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
	A = np.log2(A + A.T + 1)
	A = neighbor_ave_cpu(A, pad)
	if rp==-1:
		Q = A[:]
	else:
		Q = random_walk_cpu(A, rp)
	return [cell, Q.reshape(ngene*ngene)]

def hicluster_cpu(network, chromsize, nc, res=1000000, pad=1, rp=0.5, prct=20, ndim=20, ncpus=10):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		paras = [[cell, c, ngene, pad, rp] for cell in network]
		p = Pool(ncpus)
		result = p.map(impute_cpu, paras)
		p.close()
		index = {x[0]:j for j,x in enumerate(result)}
		Q_concat = np.array([result[index[x]][1] for x in network])
		if prct>-1:
			thres = np.percentile(Q_concat, 100 - prct, axis=1)
			Q_concat = (Q_concat > thres[:, None])
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		ndim = int(min(Q_concat.shape) * 0.2) - 1
		pca = PCA(n_components = ndim)
		R_reduce = pca.fit_transform(Q_concat)
		matrix.append(R_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce, pca.explained_variance_, pca.explained_variance_ratio_    ### add pca.explained_variance_, pca.explained_variance_ratio_

def raw_pca(network, chromsize, nc=2, res=1000000, ndim=20):
	label = network ###added by yufan
	matrix = []
	for i,c in enumerate(chromsize):
		ngene = int(chromsize[c] / res) + 1								### matrix size
		start_time = time.time()										### time start for analysis
		uptri = np.triu_indices(ngene, 1)								### get upper-triangle array without diagonal
		A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)	### iniate the matrix
		j = 0
		#print(ngene, chromsize[c], uptri)
		for cell in network:
			print(cell + '_chr' + c + '.txt')
			D = np.loadtxt(cell + '_chr' + c + '.txt')					###load sample files 
			A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()		###transfer files to matrix
			A = np.log2(A + A.T + 1)									###matrix get log2 value
			A_concat[j, :] = A[uptri]									###concatenate the matrix
			j += 1
			#print(A_concat)
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		pca = PCA(n_components = min(A_concat.shape)-1)					###PCA for cells in a certain chromosome
		A_reduce = pca.fit_transform(A_concat)
		matrix.append(A_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)						###Second PCA for all chromosomes
	#pca = PCA(n_components = 3)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, 1:(ndim + 1)])
	return kmeans.labels_, matrix_reduce

def ds_pca(network, chromsize, nc, res=1000000, ndim=20):
	matrix=[]
	for i, c in enumerate(chromsize):
		ngene = int(chromsize[c] / res)+1
		start_time = time.time()
		uptri = np.triu_indices(ngene, 1)
		A_concat = np.zeros((len(label), len(uptri[0]))).astype(float)
		j = 0
		tot = int(np.min([np.sum(np.loadtxt(cell + '_chr' + c + '.txt')[:, 2]) for cell in network]))
		for cell in network:
			D = np.loadtxt(cell + '_chr' + c + '.txt')		
			A = csr_matrix((ngene, ngene)).toarray().astype(float)
			idx = np.concatenate([[j for i in range(int(D[j,2]))] for j in range(len(D))])
			shu = np.arange(len(idx))
			np.random.shuffle(shu)
			for x in idx[shu[:tot]]:
				A[int(D[int(x),0]), int(D[int(x),1])] += 1
			A = np.log2(A + A.T + 1)
			A_concat[j, :] = A[uptri]
			j += 1
		end_time = time.time()
		print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
		pca = PCA(n_components = min(A_concat.shape)-1)
		A_reduce = pca.fit_transform(A_concat)
		matrix.append(A_reduce)
		print(c)
	matrix = np.concatenate(matrix, axis = 1)
	pca = PCA(n_components = min(matrix.shape) - 1)
	matrix_reduce = pca.fit_transform(matrix)
	kmeans = KMeans(n_clusters = nc, n_init = 200).fit(matrix_reduce[:, :ndim])
	return kmeans.labels_, matrix_reduce


