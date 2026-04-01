#Takes a genomeInformation.csv (from derep) and gtdb taxonomy tsv and merges them
#run like: python3 script genomeInformation.csv gtdbtaxonomy.tsv gtdbarchtaxonomy.tsv stats.tsv outxlsx
import sys
import pandas as pd
import os.path as path
import numpy as np
import openpyxl


#save input vars
geninfo = sys.argv[1]
gtdbtax = sys.argv[2]
gtdb_arch = sys.argv[3]
statsinfo = sys.argv[4]

#save output var
final_out = sys.argv[5]

#make pandas dfs
geninfo_df = pd.read_csv(geninfo)
gtdbtax_df = pd.read_table(gtdbtax)
statsinfo_df = pd.read_table(statsinfo)

#check for archeae
if path.exists(gtdb_arch):
   gtdbarch_df = pd.read_table(gtdb_arch)
else:
   columns = ['user_genome', 'classification', 'closest_genome_reference', 'closest_genome_reference_radius', 'closest_genome_taxonomy', 'closest_genome_ani', 'closest_genome_af', 'closest_placement_reference', 'closest_placement_radius', 'closest_placement_taxonomy', 'closest_placement_ani', 'closest_placement_af', 'pplacer_taxonomy', 'classification_method', 'note', 'other_related_references(genome_id,species_name,radius,ANI,AF)', 'msa_percent', 'translation_table', 'red_value', 'warnings']
   gtdbarch_df = pd.DataFrame(columns=columns)

#merge gtdb dataframes
gtdbtax_df = pd.concat([gtdbtax_df, gtdbarch_df], ignore_index=True, sort=False)


#add .fa ending to gtdbtax
gtdbtax_df['genome'] = gtdbtax_df['user_genome'].astype(str) + '.fa'

#make new dfs with only columns we want
geninfo_df_clean = geninfo_df.loc[:, ['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'length']]
gtdbtax_df_clean = gtdbtax_df.loc[:, ['genome', 'classification']]
statsinfo_df_clean = statsinfo_df.loc[:, ['genome', 'n_scaffolds', 'scaf_N50', 'scaf_L50', 'scaf_N90', 'scaf_L90', 'scaf_max', 'gc_avg', 'gc_std']]

#rename some of the stats columns to be less confusing
statsinfo_df_clean = statsinfo_df_clean.rename(columns={"n_scaffolds" : "scaffolds", "scaf_N50" : "L50", "scaf_L50" : "N50", "scaf_N90" : "L90", "scaf_L90" : "N90", "scaf_max" : "max_length"})

#left join on GTDB, this is a bit arbitary
joined_df = gtdbtax_df_clean.merge(geninfo_df_clean, on='genome', how='left')
joined_df = joined_df.merge(statsinfo_df_clean, on='genome', how='left')
#reorder columns
joined_df = joined_df[['genome', 'completeness', 'contamination', 'strain_heterogeneity', 'scaffolds', 'length', 'gc_avg', 'gc_std', 'N50', 'L50', 'N90', 'L90', 'max_length', 'classification']]

#save joined_df as csv
joined_df.to_excel(final_out, index=False)
#joined_df.to_csv('final_QC.csv', index=False)
