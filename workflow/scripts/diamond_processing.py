import sys
import os
import pandas as pd
import numpy as np
import plotly.express as px

#get all diamond.tsv
diamond_files = [f for f in os.listdir() if f.endswith('diamond.tsv')]

#make empty df for saving
full_df = pd.DataFrame()

#loop over the files
for file in diamond_files:

    #load in file
    diamond_df = pd.read_table(file)
    sample_name = file.replace('.norRNA.diamond.tsv', '')

    #CLEAN UP THE FILE
    #remove things after unclassified (not needed in this)
    clean_df = diamond_df.replace({'Taxonomy_string': r'unclassified.*'}, {'Taxonomy_string': ''}, regex=True)
    #Remove CO from things labeled in depth
    clean_df = clean_df.replace({'Taxonomy_string': r'cellular organisms; '}, {'Taxonomy_string': ''}, regex=True)
    #Rename CO that are only labeled to CO
    clean_df = clean_df.replace({'Taxonomy_string': r'cellular organisms'}, {'Taxonomy_string': 'Cellular Organism'}, regex=True)
    #Remove this random ncbi tax ##  (might need to change later
    clean_df = clean_df.replace({'Taxonomy_string': r'2052317'}, {'Taxonomy_string': ''}, regex=True)
    #Split tax in file
    clean_df = clean_df['Taxonomy_string'].str.split(';', expand=True)

    #change empty cells to nan
    clean_df = clean_df.replace({0: r''}, {0: np.nan}, regex=True)
    #Fill all nan with Unclassified
    clean_df = clean_df.fillna(value="Unclassified")

    #get all counts, keeping Nas for now
    cleancounts_df = clean_df[0].value_counts(dropna=False).to_frame()


#Issue starts here, sample_name can't be placed in
    #rename counts
    cleancounts_df.rename(columns={"count": "sample_name"}, inplace=True)
    #move indexes
    cleancounts_df['index_column'] = cleancounts_df.index
    #melt
    merge_df = pd.melt(cleancounts_df, id_vars=['index_column'], value_vars=['sample_name'])
    #fix name in variable column
    merge_df = merge_df.replace({'variable': r'sample_name'}, {'variable': sample_name}, regex=True)

    #add to dataframe
    full_df = pd.concat([full_df, merge_df])


#change "index_column"
#change "variable"
#change "value
full_df.rename(columns={"index_column": "Type", "variable": "Sample", "value": "Counts"}, inplace=True)

fig1 = px.histogram(full_df, x="Sample", y="Counts", color="Type", barnorm= "percent", title ="Percent of norRNA Reads Classified")
fig2 = px.histogram(full_df, x="Sample", y="Counts", color="Type", title ="Counts of norRNA Reads Classified")

fig1.write_html("norRNAclassification_percent.html")
fig2.write_html("norRNAclassification_counts.html")