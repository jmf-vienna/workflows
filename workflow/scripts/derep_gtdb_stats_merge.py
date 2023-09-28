#Takes a genomeInformation.csv (from derep) and gtdb taxonomy tsv and merges them
#run like: python3 script genomeInformation.csv gtdbtaxonomy.tsv stats.tsv outxlsx
import sys
import pandas as pd
import numpy as np
import openpyxl


#save input vars
geninfo = sys.argv[1]
gtdbtax = sys.argv[2]
statsinfo = sys.argv[3]

#save output var
final_out = sys.argv[4]

#make pandas dfs
geninfo_df = pd.read_csv(geninfo)
gtdbtax_df = pd.read_table(gtdbtax)
statsinfo_df = pd.read_table(statsinfo)

#add .fa ending to gtdbtax
gtdbtax_df['genome'] = gtdbtax_df['user_genome'].astype(str) + '.fa'

#make new dfs with only columns we want
geninfo_df_clean = geninfo_df.loc[:, ['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'length', 'N50']]
gtdbtax_df_clean = gtdbtax_df.loc[:, ['genome', 'classification']]
statsinfo_df_clean = statsinfo_df.loc[:, ['genome', 'scaf_L50', 'scaf_N90', 'scaf_L90', 'scaf_max', 'gc_avg', 'gc_std']]

#rename some of the stats columns to be less confusing
statsinfo_df_clean = statsinfo_df_clean.rename(columns={"scaf_L50" : "L50", "scaf_N90" : "N90", "scaf_L90" : "L90", "scaf_max" : "max_length"})

#left join on GTDB, this is a bit arbitary
joined_df = gtdbtax_df_clean.merge(geninfo_df_clean, on='genome', how='left')
joined_df = joined_df.merge(statsinfo_df_clean, on='genome', how='left')
#reorder columns
joined_df = joined_df[['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'length', 'gc_avg', 'gc_std', 'N50', 'L50', 'N90', 'L90', 'max_length', 'classification']]

#save joined_df as csv
joined_df.to_excel(final_out, index=False)
#joined_df.to_csv('final_QC.csv', index=False)

