import argparse
import patsy
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--columns')
parser.add_argument('--data')

args = parser.parse_args()

with open(args.columns, 'r') as file:
    columns = file.read().split()

data = patsy.demo_data(*columns)

pd.DataFrame(data).to_feather(args.data)
