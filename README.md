# MUDI
Multichannel Data Integration for single cells


MUDI is a python package developed for integration of multichannel data such as single cell Hi-C, single cell RNA-seq, single cell ATAC-seq, single cell ChIP-seq, single cell DNA methylation and so on.

The first version focuses on the integration of single cell Hi-C and single cell RNA-seq data.

Usuage:

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

###TADs

from TAD import *

#class scHiC

singleCell = scHiC()
