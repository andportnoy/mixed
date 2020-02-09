import argparse
import mixed
import numpy as np
import pandas as pd

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from sksparse.cholmod import cholesky, cholesky_AAt, analyze

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

offset = np.zeros(n)

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

cu = L.solve_L(L.apply_P(Lambdat @ ZtWy), use_LDLt_decomposition=False)
cu.tofile('cu-py.bin')

RZX = L.solve_L(L.apply_P(Lambdat @ ZtWX), use_LDLt_decomposition=False)
RZX.tofile('RZX-py.bin')

DD = XtWX - RZX.T @ RZX
DD.tofile('DD-py.bin')

# ran into an issue with using CHOLMOD here, fall back to scipy.sparse
DD = csc_matrix(DD)
b = XtWy - RZX.T @ cu
beta = spsolve(DD, b)

u = L.apply_Pt(L.solve_Lt((cu.T - RZX @ beta)[0], use_LDLt_decomposition=False))

b = Lambdat.T @ u

mu = Z.T@b + X@beta + offset

# remember to do this in sparse mode
wtres = Whalf * (y-mu)

pwrss = (wtres*wtres).sum() + (u*u).sum()

fn = float(len(mu))
ld = L.logdet()

REML = True
if REML:
    ld += cholesky(DD).logdet()
    fn -= len(beta)

deviance = ld + fn*(1. + np.log(2.*np.pi*pwrss) - np.log(fn))
np.array([deviance]).tofile('deviance-py.bin')
