import pandas as pd
from urine_dataframes import *
from composite_score import *

# =============================================================================
# LOAD & PROCESS
# =============================================================================

urine_url = 'observations/urine/urine_data.csv'
urine_raw = pd.read_csv(urine_url, index_col=0)

# Transpose: NASA files have samples as columns, we need samples as rows
urine_long      = organize_urine_df(urine_raw.transpose())
urine_processed = process_urine_df(urine_long, 'urine_alamar_656')

# =============================================================================
# PERSISTENCE TABLE
# Which markers stayed significantly changed across multiple postflight
# timepoints? Useful for writeup observations.
# =============================================================================

postflight  = urine_processed[urine_processed['timepoint'].str.startswith('R+')]

# Flag significance using the same threshold as composite scoring
postflight  = postflight.copy()
postflight['significant'] = postflight['log2fc'].abs() >= SIGNIFICANCE_THRESHOLD

significant = postflight[postflight['significant']]

persistence = (
    significant
    .groupby(['astronaut', 'marker'])
    .agg(
        times_significant = ('log2fc', 'count'),
        max_log2fc        = ('log2fc', lambda x: x.abs().max()),
        mean_log2fc       = ('log2fc', 'mean'),
        direction         = ('log2fc', lambda x: 'up' if x.mean() > 0 else 'down'),
    )
    .reset_index()
)

print("\n=== URINE MARKER PERSISTENCE (OSD-656) ===")
for astronaut, data in persistence.groupby('astronaut'):
    print(f'\n--- {astronaut} ---')
    ranked = data.sort_values(['times_significant', 'max_log2fc'], ascending=False)
    print(ranked.head(10).to_string(index=False))

# =============================================================================
# COMPOSITE DOMAIN SCORES
# =============================================================================

print("\n\n=== DOMAIN COMPOSITE SCORES ===")
timepoint_df, overall_df = calculate_composite_scores(
    urine_processed,
    DOMAIN_MAP,
    weight_map=DOMAIN_WEIGHTS,
)