import numpy as np
import pandas as pd
 
# =============================================================================
# DATA LOADING & PROCESSING
# =============================================================================
 
def organize_df(df):
    """Reshape raw NASA CSV into long format with marker, astronaut, timepoint."""
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})
    melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')
    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')
    return melted_df
 
def compute_baseline(df):
    """Average all preflight (L-) timepoints per astronaut per marker."""
    preflight = df[df['timepoint'].str.startswith('L-')]
    return (
        preflight
        .groupby(['astronaut', 'marker'])['concentration']
        .mean()
        .reset_index()
        .rename(columns={'concentration': 'baseline_mean'})
    )
 
def process_cytokine_df(df, name):
    """Compute personal baseline and log2FC for every marker. No scoring yet."""
    baseline = compute_baseline(df)
    df = df.merge(baseline, on=['astronaut', 'marker'])
 
    # Guard against log2(0)
    df = df[(df['concentration'] > 0) & (df['baseline_mean'] > 0)]
 
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    df['source'] = name
    return df
 
# =============================================================================
# SCORING
# Significance threshold: |log2FC| >= 1.0 (2x fold-change)
# Domain score at each timepoint = fraction of significant markers
# Overall domain score = mean of timepoint fractions
# Risk tier = fraction binned to 1-5
# =============================================================================
 
SIGNIFICANCE_THRESHOLD = 1.0  # |log2FC| >= 1.0 = 2x fold-change
 
def fraction_to_risk(fraction):
    """
    Convert fraction of significant markers to a 1-5 risk score.
    Thresholds represent intuitive quartile-like cutoffs:
      1 = minimal  — <10% of domain markers affected
      2 = mild     — 10-25% affected
      3 = moderate — 25-50% affected
      4 = significant — 50-75% affected
      5 = severe   — >75% of domain markers affected
    """
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                 return 5
 
def score(df, markers, domain_name=''):
    """
    Score a biological domain for each astronaut.
 
    For each postflight timepoint:
      - Count how many domain markers have |log2FC| >= threshold (significant)
      - Compute fraction = significant / total domain markers
 
    Overall score = mean fraction across all postflight timepoints.
    Risk score = fraction binned to 1-5.
 
    Also prints:
      - Per-timepoint fraction table
      - Which markers were significant and their direction
      - Overall risk scores
 
    Parameters
    ----------
    df      : processed dataframe with log2fc column
    markers : list of marker names belonging to this domain
    domain_name : label for printouts
 
    Returns
    -------
    timepoint_df  : per-astronaut × timepoint fraction and significant markers
    overall_df    : per-astronaut overall score and risk tier
    """
    # Filter to postflight only and domain markers only
    domain_df = df[
        df['timepoint'].str.startswith('R+') &
        df['marker'].isin(markers)
    ].copy()
 
    total_markers = len(markers)
 
    # Flag significant markers
    domain_df['significant'] = domain_df['log2fc'].abs() >= SIGNIFICANCE_THRESHOLD
    domain_df['direction'] = np.where(
        domain_df['log2fc'] > 0, 'up',
        np.where(domain_df['log2fc'] < 0, 'down', 'stable')
    )
 
    # Per timepoint: fraction significant + list of significant markers
    def summarize_timepoint(group):
        sig = group[group['significant']]
        return pd.Series({
            'n_significant':  len(sig),
            'n_total':        total_markers,
            'fraction':       len(sig) / total_markers,
            'up_markers':     ', '.join(sig[sig['direction'] == 'up']['marker'].tolist()),
            'down_markers':   ', '.join(sig[sig['direction'] == 'down']['marker'].tolist()),
        })
 
    timepoint_df = (
        domain_df
        .groupby(['astronaut', 'timepoint'])
        .apply(summarize_timepoint)
        .reset_index()
    )
 
    # Sort timepoints numerically
    timepoint_df['timepoint_num'] = timepoint_df['timepoint'].str.extract(r'R\+(\d+)').astype(int)
    timepoint_df = timepoint_df.sort_values(['astronaut', 'timepoint_num']).drop(columns='timepoint_num')
 
    # Overall score: mean fraction across timepoints → risk tier
    overall_df = (
        timepoint_df
        .groupby('astronaut')['fraction']
        .mean()
        .reset_index()
        .rename(columns={'fraction': 'overall_fraction'})
    )
    overall_df['risk_score'] = overall_df['overall_fraction'].apply(fraction_to_risk)
    overall_df = overall_df.sort_values('risk_score', ascending=False)
 
    # --- Printout ---
    print(f"\n{'='*60}")
    print(f"  {domain_name} (threshold: |log2FC| >= {SIGNIFICANCE_THRESHOLD})")
    print(f"{'='*60}")
 
    print("\nPer-timepoint fraction of significant markers:")
    pivot = timepoint_df.pivot_table(
        index='astronaut', columns='timepoint', values='fraction'
    )
    ordered_cols = sorted(pivot.columns, key=lambda x: int(x.replace('R+', '')))
    print(pivot[ordered_cols].round(2).to_string())
 
    print("\nSignificant markers per astronaut per timepoint:")
    for _, row in timepoint_df.iterrows():
        up   = row['up_markers']   or 'none'
        down = row['down_markers'] or 'none'
        print(f"  {row['astronaut']} {row['timepoint']}: "
              f"{int(row['n_significant'])}/{int(row['n_total'])} significant | "
              f"up: {up} | down: {down}")
 
    print("\nOverall risk scores:")
    print(overall_df.round(3).to_string(index=False))
 
    return timepoint_df, overall_df


# import numpy as np
# import pandas as pd

# # Clean dataframe for easy studies
# def organize_df(df):
#     # Modify dataframe to focus only on concentrations
#     conc = df.reset_index()
#     conc = conc[conc['index'].str.contains('concentration')]
#     conc = conc.rename(columns={'index': 'marker'})

#     # Melt dataframe into 5 different columns: marker, sample id, concentration, astronaut and timepoint
#     melted_df = conc.melt(id_vars = 'marker', var_name = 'sample_id', value_name = 'concentration')
#     melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')

#     return melted_df

# # Get the baseline mean
# def compute_baseline(df):
#     # Find whatever preflight timepoints exist in this dataset
#     preflight_mask = df['timepoint'].str.startswith('L-')
#     preflight = df[preflight_mask]
    
#     # Average them per astronaut per marker
#     baseline = preflight.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
#     baseline = baseline.rename(columns={'concentration': 'baseline_mean'})

#     return baseline












# def process_cytokine_df(df, name):
#     # Compute personal baseline from whatever preflight timepoints exist
#     baseline = compute_baseline(df)
    
#     # Merge and compute log2fc
#     df = df.merge(baseline, on=['astronaut', 'marker'])
#     df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    
#     # Score
#     df['marker_score'] = df['log2fc'].apply(score_log2fc)
    
#     # Tag which dataset this came from
#     df['source'] = name
    
#     return df

# def score(df, markers):
#     df = df[
#         df['timepoint'].str.startswith('R+') &
#         df['marker'].isin(markers)
#     ]

#     composite = (
#         df.groupby(['astronaut', 'timepoint'])['marker_score']
#         .mean()
#         .reset_index()
#     )
#     composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
#     composite = composite.sort_values(['astronaut', 'timepoint_num'])

#     # Print per-timepoint scores per astronaut
#     print("\nPer-timepoint scores:")
#     per_tp = (
#         composite
#         .pivot_table(index='astronaut', columns='timepoint', values='marker_score')
#     )
#     ordered_cols = sorted(per_tp.columns, key=lambda x: int(x.replace('R+','')))
#     per_tp = per_tp[ordered_cols].round(2)
#     print(per_tp.to_string())

#     overall_scores = (
#         composite
#         .groupby('astronaut')['marker_score']
#         .mean()
#         .reset_index()
#         .rename(columns={'marker_score': 'overall_score'})
#         .sort_values('overall_score', ascending=False)
#     )

#     # Print overall weighted scores
#     print("\nOverall weighted score:")
#     print(overall_scores.round(2).to_string(index=False))

#     return composite, overall_scores
