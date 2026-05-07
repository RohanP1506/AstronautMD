import pandas as pd

'''Making immune_eve dataframe'''

immune_eve = pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv', index_col=0).transpose()

# Cleaning immune_eve data to only focus on concentrations
immune_eve = immune_eve.reset_index()
immune_eve = immune_eve.rename(columns={'index': 'Sample Name'})
immune_eve_conc = immune_eve[immune_eve["Sample Name"].str.contains('concentration')]

'''Organizing dataframe for easier calculations'''

# Step 1 — melt
df_long = immune_eve_conc.melt(
    id_vars='Sample Name',      # keep this column as-is
    var_name='sample_id',       # all the column headers become values in this column
    value_name='concentration'  # all the numbers go here
)

# Step 2 — split C001_serum_L-3 into separate columns
df_long[['astronaut', 'fluid', 'timepoint']] = df_long['sample_id'].str.split('_', n=2, expand=True)

# Step 3 — rename for clarity
df_long = df_long.rename(columns={'Sample Name': 'marker'})

'''Compute personal baselines to compare with log2'''

import numpy as np

# Filter to pre-flight timepoints only
baseline_tps = ['L-3', 'L-44', 'L-92']
df_baseline = df_long[df_long['timepoint'].isin(baseline_tps)]

# Average the 3 pre-flight timepoints per astronaut per marker
baseline = df_baseline.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
baseline = baseline.rename(columns={'concentration': 'baseline_mean'})

# Attach each astronaut's baseline to every row
df_long = df_long.merge(baseline, on=['astronaut', 'marker'])

# Compute log2 fold change vs personal baseline
df_long['log2fc'] = np.log2(df_long['concentration'] / df_long['baseline_mean'])

def score_log2fc(x):
    x = abs(x)
    if x < 0.5:   return 0   # stable
    elif x < 1.0: return 1   # mild shift
    elif x < 1.5: return 2   # moderate shift
    else:         return 3   # large shift

df_long['marker_score'] = df_long['log2fc'].apply(score_log2fc)

print(df_long[['marker', 'astronaut', 'timepoint', 'concentration', 'log2fc', 'marker_score']].head(20))

# Next step is persistence weighting #
