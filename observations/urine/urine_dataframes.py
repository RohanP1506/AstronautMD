import numpy as np
import pandas as pd

# =============================================================================
# DATA LOADING & PROCESSING
# =============================================================================

def organize_urine_df(df):
    """
    Reshape a transposed NASA urine CSV into long format.
    Only keeps concentration columns. Filters out zero and NaN concentrations.
    Output columns: marker, sample_id, concentration, astronaut, timepoint
    """
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})

    melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')

    # Filter out zeros and NaNs before any math
    melted_df['concentration'] = pd.to_numeric(melted_df['concentration'], errors='coerce')
    melted_df = melted_df[melted_df['concentration'] > 0].dropna(subset=['concentration'])

    # Parse astronaut ID and timepoint
    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_urine_(.*)')
    melted_df = melted_df.dropna(subset=['astronaut', 'timepoint'])

    return melted_df.reset_index(drop=True)


def compute_baseline(df):
    """
    Compute each astronaut's personal baseline per marker as the mean
    concentration across all preflight timepoints (L-).
    Excludes zero baselines to avoid division by zero in log2FC.
    """
    preflight = df[df['timepoint'].str.startswith('L-')]
    baseline = (
        preflight
        .groupby(['astronaut', 'marker'])['concentration']
        .mean()
        .reset_index()
        .rename(columns={'concentration': 'baseline_mean'})
    )
    # Drop zero baselines — can't compute log2FC against zero
    baseline = baseline[baseline['baseline_mean'] > 0]
    return baseline


def process_urine_df(df, name):
    """
    Merge baseline into the dataframe and compute log2FC per marker per
    timepoint relative to that astronaut's personal preflight baseline.
    No scoring here — scoring happens in score().
    """
    baseline = compute_baseline(df)

    # Inner merge drops any astronaut/marker combos without a valid baseline
    df = df.merge(baseline, on=['astronaut', 'marker'])

    # Both concentration and baseline_mean are guaranteed > 0 from above filters
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    df['source'] = name

    return df.reset_index(drop=True)