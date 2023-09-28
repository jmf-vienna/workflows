#Take genome coverage generated from LISC outputs and pivot them into a table that could be merged with other tables
##run like: python3 script coveragestats.csv outxlsx

import sys
import pandas as pd
import openpyxl

#save input name
covstats = sys.argv[0]
final_out = sys.argv[1]

#open table as dataframe
covstats_df = pd.read_csv(covstats)

#reduce table to only columns of interest
covstats_df_clean = covstats_df.loc[:, ['genome', 'library', 'average_cov', 'percent_bases_mapped']]

#shorten some column names
covstats_df_clean = covstats_df_clean.rename(columns={"average_cov" : "avg_cov", "percent_bases_mapped" : "perc_mapped"})

#mop up the library names
covstats_df_clean['library'] = covstats_df_clean['library'].str.replace('_pass.interleaved.fastq.gz', '')

#pivot table
covstats_piv = covstats_df_clean.pivot(index='genome', columns='library', values='avg_cov')


#print out
covstats_piv.to_excel(final_out, index=True)