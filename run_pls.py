import argparse
import mixed
import numpy as np
import pandas as pd

from scipy.sparse import csc_matrix
from sksparse.cholmod import cholesky_AAt

parser = argparse.ArgumentParser()
parser.add_argument('--data')
parser.add_argument('--formula')
parser.add_argument('--randomdata')

args = parser.parse_args()

data = pd.read_feather(args.data)
with open(args.formula, 'r') as file:
    formula = file.read().strip()

X, Z, Lambdat, y, thfun = mixed.get_matrices(data, formula)

n = len(y)

Whalf = np.eye(n)
Whalf.tofile('Whalf-py.bin')

WX = Whalf @ X
WX.tofile('WX-py.bin')

Wy = Whalf @ y
Wy.tofile('Wy-py.bin')

ZtW = Z @ Whalf
ZtW.tofile('ZtW-py.bin')

XtWX = X.T @ WX
XtWX.tofile('XtWX-py.bin')

XtWy = WX.T @ Wy
XtWy.tofile('XtWy-py.bin')

ZtWX = ZtW @ WX
ZtWX.tofile('ZtWX-py.bin')

ZtWy = ZtW @ Wy
ZtWy.tofile('ZtWy-py.bin')

DD = XtWX

LambdatZtW = csc_matrix(Lambdat @ ZtW)

L = cholesky_AAt(LambdatZtW, beta=1)
L.L().todense().tofile('L-py.bin')

newtheta = np.fromfile(args.randomdata)[:len(Lambdat.data)]

# deviance function calculations
Lambdat.data[:] = thfun(newtheta)
Lambdat.todense().tofile('Lambdat-new-py.bin')

LambdatZtW = csc_matrix(Lambdat @ ZtW)
L = cholesky_AAt(LambdatZtW, beta=1)
L.L().todense().tofile('L-new-py.bin')
