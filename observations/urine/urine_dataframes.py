import numpy as np
import pandas as pd

# 1. Adjusted organization function for Urine naming conventions
def organize_urine_df(df):
    # Modify dataframe to focus only on concentrations
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})

    # Melt dataframe: marker stays as ID, sample names become a column
    melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')
    
    # Updated Regex: Looks for "_urine_" instead of "_serum_"
    # Pattern: (Astronaut ID)_urine_(Timepoint)
    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_urine_(.*)')

    return melted_df

def score_log2fc(x):
    x = abs(x)
    if x < 0.5: return 0
    elif x < 1.0: return 1
    elif x < 1.5: return 2
    else: return 3

def compute_baseline(df):
    # Standard I4 preflight timepoints start with 'L-' (e.g., L-92, L-44, L-3)
    preflight_mask = df['timepoint'].str.startswith('L-')
    preflight = df[preflight_mask]
    
    baseline = preflight.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
    baseline = baseline.rename(columns={'concentration': 'baseline_mean'})
    return baseline

def process_cytokine_df(df, name):
    baseline = compute_baseline(df)
    df = df.merge(baseline, on=['astronaut', 'marker'])
    
    # Calculate Log2 Fold Change relative to the astronaut's own baseline
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    df['marker_score'] = df['log2fc'].apply(score_log2fc)
    df['source'] = name
    return df