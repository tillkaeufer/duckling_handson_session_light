from dust_extinction.parameter_averages import G23 as ext_model
import glob
import numpy as np
import time
import os
from scipy import interpolate
from scipy.optimize import nnls
import json
import uuid
import multiprocessing as mp

from PyAstronomy import pyasl
import corner

import ultranest
import matplotlib.pyplot as plt


import sys
import importlib

from spectres.spectral_resampling_numba import spectres_numba  as spectres

from utils import *


from matplotlib.lines import Line2D
import pickle 
import matplotlib as mpl
import matplotlib.pyplot as plt
from ast import literal_eval

import argparse

print('All packages are loaded')
print('Congrats this means that everthing is installed correctly!')
