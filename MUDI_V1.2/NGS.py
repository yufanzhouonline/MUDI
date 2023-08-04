import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb

from abc import ABC, abstractmethod

class NGS(ABC):
    '''
    NGS class is the Abstract Base Class used for Next Generation Sequencing data.
    The sub-class could be RNAseq, ATACseq, ChIPseq or HiC and so on.
    '''
    
    @abstractmethod
    def load(self):
        pass


