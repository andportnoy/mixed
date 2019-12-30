import argparse
import mixed
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--data')
parser.add_argument('--formula')

args = parser.parse_args()

data = pd.read_feather(args.data)
with open(args.formula, 'r') as file:
    formula = file.read().strip()

X, Z, Lambdat, y = mixed.get_matrices(data, formula)

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

LambdatZtW = Lambdat @ ZtW
I = np.eye(Lambdat.shape[0])

L = np.linalg.cholesky(LambdatZtW @ LambdatZtW.T + I)
L.tofile('L-py.bin')
