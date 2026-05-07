import numpy as np
import pandas as pd
from cytokine_dataframe import *

# Dataset straight from open-source
immune_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv', index_col=0).transpose())

# Clean up dataset
immune_eve_long = process_cytokine_df(immune_eve, 'immune_eve')

# Make observations on significant changes in concentration before and after flight
postflight = immune_eve_long[immune_eve_long['timepoint'].str.startswith('R+')]
significant = postflight[postflight['marker_score'] >= 1]

# Make a dataframe of the persistence of significant changes in certain markers over timepoints
persistence = significant.groupby(['astronaut', 'marker']).agg(
    times_significant = ('marker_score', 'count'), # count rows
    max_score = ('marker_score', 'max'), # highest score
    mean_log2fc = ('log2fc', 'mean') # average log2fc
).reset_index()

# # Rank persistency between all astronauts
# for astronaut, data in persistence.groupby('astronaut'):
#     print(f'\n=== {astronaut} ===')
#     ranked = data.sort_values('times_significant', ascending=False)
#     print(ranked)
