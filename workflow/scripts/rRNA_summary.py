import sys
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots


def split_name(name):
    return pd.Series(name.split(".", 1))

All_stats = pd.read_table("stats.tsv")

All_stats_clean = All_stats.replace({'file': r'^.*/'}, {'file': ''}, regex=True)
All_stats_clean[['Sample', 'Type']] = All_stats_clean['file'].apply(split_name)

All_stats_clean = All_stats_clean.replace({'Type': r'rRNA.fq.gz'}, {'Type': ''}, regex=True)
All_stats_clean = All_stats_clean.replace({'Type': r'no'}, {'Type': 'non-rRNA'}, regex=True)
All_stats_clean = All_stats_clean.replace({'Type': r'cleaned.fastq.gz'}, {'Type': 'full sample'}, regex=True)
All_stats_clean = All_stats_clean.replace({'Type': r'rRNA.800.fq.gz'}, {'Type': '_800'}, regex=True)


fig1 = px.bar(All_stats_clean, x="Sample", y="num_seqs", color='Type', barmode='group')


id_stats = All_stats_clean[All_stats_clean["Type"].str.contains("800|full") == False]

fig2 = px.histogram(id_stats, x="Sample", y="num_seqs", color="Type", barnorm= "percent", title ="Percent of Ribosomal Sequences")
fig3 = px.histogram(id_stats, x="Sample", y="num_seqs", color="Type", title ="Number of Ribosomal Sequences")

fig1.write_html("All_stats.html")
fig2.write_html("Sequences_percent.html")
fig3.write_html("Sequences_counts.html")

