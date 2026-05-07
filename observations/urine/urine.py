import numpy as np
import pandas as pd
from urine_dataframes import *


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