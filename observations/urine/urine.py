import numpy as np
import pandas as pd
from urine_dataframes import *
from composite_score import *
from urine_viz import generate_all_figures


# 2. Load the OSD-656 Urine Dataset
urine_url = 'observations/urine/urine_data.csv'
urine_raw = pd.read_csv(urine_url, index_col=0)

# 3. Run Pipeline
# We transpose because NASA 'TRANSFORMED' files usually have samples as rows, 
# and your logic expects samples as columns before the melt.
urine_long = organize_urine_df(urine_raw.transpose())
urine_processed = process_cytokine_df(urine_long, 'urine_alamar_656')

# 4. Analyze Persistence (R+ timepoints)
postflight = urine_processed[urine_processed['timepoint'].str.startswith('R+')]
significant = postflight[postflight['marker_score'].abs() >= 2]

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

# 4. Generate Composite Scores
composite_scores = calculate_composite_scores(urine_processed, DOMAIN_MAP)

# 5. Look at the specific domains for all astronauts
print("\n=== DOMAIN COMPOSITE SCORES (0-5 SCALE) ===")
# Pivot for easier viewing: Columns are domains, rows are astronaut/timepoint
final_view = composite_scores.pivot_table(
    index=['astronaut', 'timepoint'], 
    columns='domain', 
    values='composite_score'
).reset_index()

print(final_view.sort_values(['astronaut', 'timepoint']))

# Check how many post-flight (R+) timepoints had a composite score > 2.0
postflight_domains = composite_scores[composite_scores['timepoint'].str.startswith('R+')]
domain_persistence = postflight_domains[postflight_domains['composite_score'] >= 2.0]

print("\n--- Persistent Domain Elevations (R+ Timepoints) ---")
print(domain_persistence)


generate_all_figures(urine_processed, composite_scores, output_dir='observations/urine/figures')