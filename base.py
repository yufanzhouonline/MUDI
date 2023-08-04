import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
import json

from version import *

WORK_DIR = os.getcwd()

def set_dir(work_dir = WORK_DIR):
    os.chdir(work_dir)
    print('Set current working directory as: ' + work_dir)

def get_dir():
    work_dir = os.getcwd()
    print('Current working directory is: ' + work_dir)

class MUDI(object):
    '''
    set version and parameters
    '''
    __version__ = '2.1'
    
    def __init__(self, config_path = 'config.json'):
        '''
        read json configuration file
        config: store configure parameters
        '''
        with open(config_path) as f:
            self.config = json.load(f)
    

class Base(object):
    '''
    Class Base is the Base of all classes, changing of this class will lead to changes in other classes.
    The sub-class could be any other clase based on this class such as Data, Analysis, Plot.
    '''
    def __init__(self):
        '''
        initiate the working environment
        '''
        pass


