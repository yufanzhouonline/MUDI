from base import *

COLOR_MAP = ['#13EAC988', '#FF7F5088', '#00640088', '#ED0DD988', '#DAA52088', '#65370088', '#C1F80A88', '#88888888', '#E6DAA688']
class Plot(Base):
    '''
    Class Plot is the class for plotting figures.
    The sub-class could be some specific figure class.
    '''
    
    def __init__(self):
        pass
    
    @staticmethod
    def all2D(data, cell, pc1=0, pc2=1, loc='upper right'):
        '''
        plot 2D for all single cells
        '''
        plt.figure()
        if cell.shape[0] == 1:
            p1, = plt.plot(data[1][0:cell[0], pc1],data[1][0:cell[0], pc2], '.', color=COLOR_MAP[0])
            plt.legend([p1], list(cell.index), loc=loc)
        elif cell.shape[0] == 2:
            p1, = plt.plot(data[1][0:cell[0], pc1],data[1][0:cell[0], pc2], '.', color=COLOR_MAP[0])
            p2, = plt.plot(data[1][cell[0]:(cell[0] + cell[1]), pc1],data[1][cell[0]:(cell[0] + cell[1]), pc2], '.', color=COLOR_MAP[1])
            plt.legend([p1, p2], list(cell.index), loc=loc)
        elif cell.shape[0] == 3:
            p1, = plt.plot(data[1][0:cell[0], pc1],data[1][0:cell[0], pc2], '.', color=COLOR_MAP[0])
            p2, = plt.plot(data[1][cell[0]:(cell[0] + cell[1]), pc1],data[1][cell[0]:(cell[0] + cell[1]), pc2], '.', color=COLOR_MAP[1])
            p3, = plt.plot(data[1][(cell[0] + cell[1]):(cell[0] + cell[1] + cell[2]), pc1],data[1][(cell[0] + cell[1]):(cell[0] + cell[1] + cell[2]), pc2], '.', color=COLOR_MAP[2])
            plt.legend([p1, p2, p3], list(cell.index), loc=loc)
        plt.xlabel('PC' + str(pc1 + 1), fontsize=18)
        plt.ylabel('PC' + str(pc2 + 1), fontsize=18)
        #plt.axis('square')
        plt.show()
    
    @staticmethod
    def all3D(data, cell, pc1=0, pc2=1, pc3=2, loc='upper left'):
        '''
        plot 3D cluster for all single cells
        '''
        plt.figure()
        ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
        ####separate the data point to three parts, marked them with different color
        if cell.shape[0] == 1:
            p1 = ax.scatter(data[1][:, pc1], data[1][:, pc2], data[1][:, pc3], color=COLOR_MAP[0])  ###draw the data points
            plt.legend([p1], list(cell.index), loc=loc, fontsize=8)
        elif cell.shape[0] == 2:
            p1 = ax.scatter(data[1][0:cell[0], pc1], data[1][0:cell[0], pc2], data[1][0:cell[0], pc3], color=COLOR_MAP[0])  ###draw the data points
            p2 = ax.scatter(data[1][cell[0]:(cell[0] + cell[1]), pc1], data[1][cell[0]:(cell[0] + cell[1]), pc2], data[1][cell[0]:(cell[0] + cell[1]), pc3], color=COLOR_MAP[1])
            plt.legend([p1, p2], list(cell.index), loc=loc, fontsize=8)
        elif cell.shape[0] == 3:
            p1 = ax.scatter(data[1][0:cell[0], pc1], data[1][0:cell[0], pc2], data[1][0:cell[0], pc3], color=COLOR_MAP[0])  ###draw the data points
            p2 = ax.scatter(data[1][cell[0]:(cell[0] + cell[1]), pc1], data[1][cell[0]:(cell[0] + cell[1]), pc2], data[1][cell[0]:(cell[0] + cell[1]), pc3], color=COLOR_MAP[1])
            p3 = ax.scatter(data[1][(cell[0] + cell[1]):(cell[0] + cell[1] + cell[2]), pc1], data[1][(cell[0] + cell[1]):(cell[0] + cell[1] + cell[2]), pc2], data[1][(cell[0] + cell[1]):(cell[0] + cell[1] + cell[2]), pc3], color=COLOR_MAP[1])
            plt.legend([p1, p2, p3], list(cell.index), loc=loc, fontsize=8)
        ax.set_xlabel('PC' + str(pc1 + 1))
        ax.set_ylabel('PC' + str(pc2 + 1))
        ax.set_zlabel('PC' + str(pc3 + 1))
        ax.view_init(18, 161)
        plt.show()
        
    @staticmethod
    def cluster_number(X):
        '''
        determine the optimal cluster number by Silhouette Coefficient
        '''
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
        
    @staticmethod
    def cluster2D(X, yhat, pc1=0, pc2=1):
        '''
        plot cluster 2D figure
        if pc1 = 0 and pc2 = 1, then plot PC1 and PC2
        if pc1 = 1 and pc2 = 2, then plot PC2 and PC3
        if pc1 = 0 and pc2 = 2, then plot PC1 and PC3
        '''
        plt.figure()
        cluster_array = np.unique(yhat)
        for each_cluster in cluster_array:
            row_ix = np.where(yhat == each_cluster)
            plt.scatter(X[row_ix, pc1], X[row_ix, pc2])
        plt.xlabel('PC' + str(pc1 + 1), fontsize=18)
        plt.ylabel('PC' + str(pc2 + 1), fontsize=18)
        plt.show()
        
    @staticmethod
    def cluster3D(X, yhat, pc1=0, pc2=1, pc3=2):
        '''
        plot cluster 3D figure
        '''
        plt.figure()
        cluster_array = np.unique(yhat)
        ax = plt.subplot(111, projection='3d')  ###create a 3D plot project
        for each_cluster in cluster_array:
            row_ix = np.where(yhat == each_cluster)
            ax.scatter(X[row_ix, 0], X[row_ix, 1], X[row_ix, 2])
        ax.set_xlabel('PC' + str(pc1 + 1))
        ax.set_ylabel('PC' + str(pc2 + 1))
        ax.set_zlabel('PC' + str(pc3 + 1))
        ax.view_init(18, 161)
        plt.show()



