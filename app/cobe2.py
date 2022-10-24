import numpy as np
from scipy.optimize import least_squares
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

data = pd.read_csv("src/cobe/cobeplot.csv", header=None, delim_whitespace=True)
vbi = np.asarray(data[0])
Ii = np.asarray(data[1])
vi = 30 * vbi

def r(b, vi, Ii):
    return ((0.0014745*pow(vi,3))/(np.exp(b*vi)-1) - Ii)**2

b = least_squares(r, 0.048, args=(vi, Ii))
print(b.x)
