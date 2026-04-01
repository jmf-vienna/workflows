#Takes a genomeInformation.csv (from derep) and gtdb taxonomy tsv and merges them
#run like: python3 script stats_file coverage_file outxlsx
import sys
import pandas as pd
import numpy as np
import openpyxl


#save input vars
statsinfo = sys.argv[1]
covinfo = sys.argv[2]

#save output var
final_out = sys.argv[3]

#make pandas dfs
statsinfo_df = pd.read_excel(statsinfo)
covinfo_df = pd.read_csv(covinfo)


#Fix covinfo file
#clean up
cov_df_clean = covinfo_df.rename(columns={"average_cov" : "avg_cov", "percent_bases_mapped" : "perc_mapped"})
cov_df_clean['library'] = cov_df_clean['library'].str.replace('.interleave.fastq.gz', '')
#remove total_bases and mapped_bases
cov_df_clean = cov_df_clean[['genome','library', 'avg_cov', 'perc_mapped']]
#pivot
pivcov_df = cov_df_clean.pivot(index="genome", columns='library')

#fix column names
pivcov_df.columns = ['_'.join(str(s).strip() for s in col if s) for col in pivcov_df.columns]


#merge allstats with cov

joined_df = statsinfo_df.merge(pivcov_df, on='genome', how='left')

#save joined_df as csv
joined_df.to_excel(final_out, index=False)

