import argparse
import mixed
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--data')
parser.add_argument('--formula')
parser.add_argument('--X')
parser.add_argument('--Z')
parser.add_argument('--Lambdatx')
parser.add_argument('--Lambdati')
parser.add_argument('--Lambdatp')

args = parser.parse_args()

data = pd.read_feather(args.data)
with open(args.formula, 'r') as file:
    formula = file.read().strip()

X, Z, Lambdat, _, _ = mixed.get_matrices(data, formula)

X.tofile(args.X)
Z.todense().tofile(args.Z)
Lambdat.data.tofile(args.Lambdatx)
Lambdat.indices.tofile(args.Lambdati)
Lambdat.indptr.tofile(args.Lambdatp)
