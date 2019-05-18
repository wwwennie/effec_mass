#!/usr/bin/env python3

# common packages and constants 

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

#TeX formatting in plots
from matplotlib import rc
rc('text', usetex=True)

# in SI units
from scipy.constants import c, epsilon_0, m_e, Rydberg, Boltzmann, e, hbar

bohr2ang = 0.529177249
ryd2ev = 13.605698066
ryd2har = 0.5
ev2cm1 = 8065.54429
ev2hz = 241.79893*1e12
ev2joule = 1.60218e-19
ang2m = 1e-10
