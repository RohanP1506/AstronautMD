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

def score_log2fc(x): # Scores on a 5 point scale
    x = abs(x) # treat increases and decreases equally
    if x < 0.5:
        return 1 # stable - no meaningful change
    elif x < 1.0:
        return 2 # mild shift - less than 2x change
    elif x < 1.5:
        return 3 # moderate shift - 2x to ~3x change
    elif x < 2.0:
        return 4 # large shift - more than 3x change
    else:
        return 5 # severe - greater than 4x change

# Get the baseline mean
def compute_baseline(df):
    # Find whatever preflight timepoints exist in this dataset
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

def score(df, markers):
    df = df[
        df['timepoint'].str.startswith('R+') &
        df['marker'].isin(markers)
    ]

    composite = (
        df.groupby(['astronaut', 'timepoint'])['marker_score']
        .mean()
        .reset_index()
    )
    composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
    composite = composite.sort_values(['astronaut', 'timepoint_num'])

    overall_scores = (
        composite
        .groupby('astronaut')['marker_score']
        .mean()
        .reset_index()
        .rename(columns={'marker_score': 'overall_score'})
        .sort_values('overall_score', ascending=False)
    )

    print(overall_scores)
    return composite, overall_scores
