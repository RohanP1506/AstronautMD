import numpy as np
import pandas as pd

# Clean dataframe for easy studies
def organize_df(df):
    # Modify dataframe to focus only on concentrations
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})

    # Melt dataframe into 5 different columns: marker, sample id, concentration, astronaut and timepoint
    melted_df = conc.melt(id_vars = 'marker', var_name = 'sample_id', value_name = 'concentration')
    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')

    return melted_df

def score_log2fc(x):
    x = abs(x)        # treat increases and decreases equally
    if x < 0.5:
        return 0      # stable — no meaningful change
    elif x < 1.0:
        return 1      # mild shift — less than 2x change
    elif x < 1.5:
        return 2      # moderate shift — 2x to ~3x change
    else:
        return 3      # large shift — more than 3x change

# Get the baseline mean
def compute_baseline(df):
    # Automatically find whatever preflight timepoints exist in this dataset
    preflight_mask = df['timepoint'].str.startswith('L-')
    preflight = df[preflight_mask]
    
    # Average them per astronaut per marker
    baseline = preflight.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
    baseline = baseline.rename(columns={'concentration': 'baseline_mean'})
    
    return baseline

def process_cytokine_df(df, name):
    # Compute personal baseline from whatever preflight timepoints exist
    baseline = compute_baseline(df)
    
    # Merge and compute log2fc
    df = df.merge(baseline, on=['astronaut', 'marker'])
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    
    # Score
    df['marker_score'] = df['log2fc'].apply(score_log2fc)
    
    # Tag which dataset this came from
    df['source'] = name
    
    return df

immune_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv', index_col=0).transpose())
immune_alamar = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum.immune.AlamarPanel_TRANSFORMED.csv', index_col=0).transpose())
cardio_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_cardiovascular_EvePanel_TRANSFORMED.csv', index_col=0).transpose())

immune_eve_long = process_cytokine_df(immune_eve, 'immune_eve')
immune_alamar_long = process_cytokine_df(immune_alamar, 'immune_alamar')
cardio_eve_long = process_cytokine_df(cardio_eve, 'cardio_eve')

# Make observations on significant changes in concentration before and during/after flight
postflight = immune_eve_long[immune_eve_long['timepoint'].str.startswith('R+')]
significant = postflight[postflight['marker_score'] >= 1]

# Group datasets by astronaut
astro1 = significant[significant['astronaut'] == 'C001']
astro2 = significant[significant['astronaut'] == 'C002']
astro3 = significant[significant['astronaut'] == 'C003']
astro4 = significant[significant['astronaut'] == 'C004']
astronauts = [astro1, astro2, astro3, astro4]

# go through list of astronauts
# go through each marker and count how many times its significant through each timepoint after flight
    # If once then not significant, twice maybe, thrice probably

