#Takes a genomeInformation.csv (from derep) and gtdb taxonomy tsv and merges them
#run like: python3 script genomeInformation.csv gtdbtaxonomy.tsv
import sys
import pandas as pd
import numpy as np
import openpyxl


#save input vars
geninfo = sys.argv[1]
gtdbtax = sys.argv[2]

#save output var
final_out = sys.argv[3]

#make pandas dfs
geninfo_df = pd.read_csv(geninfo)
gtdbtax_df = pd.read_table(gtdbtax)


#add .fa ending to gtdbtax
gtdbtax_df['genome'] = gtdbtax_df['user_genome'].astype(str) + '.fa'

#make new dfs with only columns we want
geninfo_df_clean = geninfo_df.loc[:, ['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'length', 'N50']]
gtdbtax_df_clean = gtdbtax_df.loc[:, ['genome', 'classification']]


#left join on GTDB, this is a bit arbitary
joined_df = gtdbtax_df_clean.merge(geninfo_df_clean, on='genome', how='left')

#reorder columns
joined_df = joined_df[['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'length', 'N50', 'classification']]

#save joined_df as csv
joined_df.to_excel(final_out, index=False)
#joined_df.to_csv('final_QC.csv', index=False)

