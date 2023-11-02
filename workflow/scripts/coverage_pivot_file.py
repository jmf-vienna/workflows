import sys
import pandas as pd
import openpyxl

infile = sys.argv[1]
outfile = sys.argv[2]

cov_df = pd.read_csv(infile)

pivcov_df = cov_df.pivot(index="genome", columns='library')

pivcov_df.columns = ['_'.join(str(s).strip() for s in col if s) for col in pivcov_df.columns]

pivcov_df.to_csv(outfile)
