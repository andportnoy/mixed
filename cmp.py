import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-a')
parser.add_argument('-b')
args = parser.parse_args()

a = np.fromfile(args.a, dtype=np.float64)
b = np.fromfile(args.b, dtype=np.float64)

assert np.array_equal(a, b)
