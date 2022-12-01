###Input genome
from Genome import *

###for batch cells
from HiC import *
from ATACseq import *
from ChIPseq import *
from RNAseq import *

###for single cells
from scHiC import *
from scATACseq import *
from scChIPseq import *
from scRNAseq import *

###others
from TAD import *

cell = scHiC()
cell.cellLabel(labelPath = './data/samples/cellLabel.txt')
cell.cellCluster(clusterPath = './data/samples/clusterList.txt')
