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

# 2. Load the OSD-656 Urine Dataset
urine_url = 'https://osdr.nasa.gov/geode-py/ws/studies/OSD-656/download?source=datamanager&file=LSDS-64_Multiplex_urine.immune.AlamarPanel_TRANSFORMED.csv'
urine_raw = pd.read_csv(urine_url, index_col=0)

# 3. Run Pipeline
# We transpose because NASA 'TRANSFORMED' files usually have samples as rows, 
# and your logic expects samples as columns before the melt.
urine_long = organize_urine_df(urine_raw.transpose())
urine_processed = process_cytokine_df(urine_long, 'urine_alamar_656')

# 4. Analyze Persistence (R+ timepoints)
postflight = urine_processed[urine_processed['timepoint'].str.startswith('R+')]
significant = postflight[postflight['marker_score'] >= 1]

persistence = significant.groupby(['astronaut', 'marker']).agg(
    times_significant = ('marker_score', 'count'),
    max_score = ('marker_score', 'max'),
    mean_log2fc = ('log2fc', 'mean')
).reset_index()

# 5. Output Rankings
print("\n=== URINE IMMUNE PERSISTENCE (OSD-656) ===")
for astronaut, data in persistence.groupby('astronaut'):
    print(f'\n--- {astronaut} ---')
    # Rank by how many post-flight timepoints remained significant
    ranked = data.sort_values(['times_significant', 'max_score'], ascending=False)
    print(ranked.head(10))