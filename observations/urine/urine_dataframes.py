import numpy as np

def organize_urine_df(df):
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})

    melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')
    
    # --- ADDED: Filter out 0 and NA values immediately ---
    # This ensures no 0s or NaNs propagate to math functions
    melted_df = melted_df[melted_df['concentration'] > 0].dropna(subset=['concentration'])

    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_urine_(.*)')

    return melted_df

def score_log2fc(x):
    if x > 2.0:   return  5   # severe upregulation
    elif x > 1.5: return  4
    elif x > 1.0: return  3
    elif x > 0.5: return  2
    elif x > 0.25: return 1
    elif x >= -0.25: return 0 # stable
    elif x >= -0.5: return -1
    elif x >= -1.0: return -2
    elif x >= -1.5: return -3
    elif x >= -2.0: return -4
    else:           return -5  # severe downregulation

def compute_baseline(df):
    preflight_mask = df['timepoint'].str.startswith('L-')
    preflight = df[preflight_mask]
    
    # groupby().mean() automatically ignores NaNs, 
    # but since we filtered 0s in organize_urine_df, the baseline will be clean.
    baseline = preflight.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
    baseline = baseline.rename(columns={'concentration': 'baseline_mean'})
    
    # --- ADDED: Ensure baseline_mean is not 0 (to avoid division by zero) ---
    baseline = baseline[baseline['baseline_mean'] > 0]
    
    return baseline

def process_cytokine_df(df, name):
    baseline = compute_baseline(df)
    
    # We use an 'inner' merge by default so that if an astronaut/marker 
    # combo didn't have a valid baseline, it's dropped from the final analysis.
    df = df.merge(baseline, on=['astronaut', 'marker'])
    
    # Calculate Log2 Fold Change
    # Because we filtered concentration > 0 and baseline_mean > 0, 
    # this will no longer produce -inf or inf.
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    
    df['marker_score'] = df['log2fc'].apply(score_log2fc)
    df['source'] = name
    return df