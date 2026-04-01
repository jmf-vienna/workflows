#Takes a final_QC and coverM outputs and merges them
#run like: python3 script final_QC.xlsx coverM_mean_output.tsv coverM_relative_output.tsv final_QC_higQ.xlsx
import sys
import pandas as pd
import numpy as np
import openpyxl


#save input vars
stats = sys.argv[1]
coverm_mean = sys.argv[2]
coverm_abund = sys.argv[3]

#save output var
final_out = sys.argv[4]

#make pandas dfs
stats_df = pd.read_excel(stats)
coverm_mean_df = pd.read_table(coverm_mean)
coverm_abund_df = pd.read_table(coverm_abund)

#merge coverMs
coverm_df = coverm_abund_df.merge(coverm_mean_df, on='Genome', how='left')

#rename "Genome" to "genome
coverm_df = coverm_df.rename(columns={"Genome": "genome"})

#remove .fa from genomes
stats_df['genome'] = stats_df['genome'].str.replace(r'.fa', '')

#merge all using stats as reference
final_df = stats_df.merge(coverm_df, on='genome', how='outer')

#remove interleave.fastq.gz from all columns
final_df.columns=final_df.columns.str.replace('.interleave.fastq.gz','')

#save joined_df as csv
final_df.to_excel(final_out, index=False)
