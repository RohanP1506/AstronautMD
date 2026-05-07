from cytokine_dataframe import *

# Dataset straight from open-source
immune_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv', index_col=0).transpose())

# Clean up dataset
immune_eve_long = process_cytokine_df(immune_eve, 'immune_eve')

# Make observations on significant changes in concentration before and after flight
postflight = immune_eve_long[immune_eve_long['timepoint'].str.startswith('R+')]
significant = postflight[postflight['marker_score'] >= 1]

# Make a dataframe of the persistence of significant changes in certain markers over timepoints
persistence = significant.groupby(['astronaut', 'marker']).agg(
    times_significant = ('marker_score', 'count'), # count rows
    max_score = ('marker_score', 'max'), # highest score
    mean_log2fc = ('log2fc', 'mean') # average log2fc
).reset_index()

### BEGIN CLASSIFICATION ###

immune_regulation_markers = [
    'il_2_concentration_picogram_per_milliliter',    # T-cell proliferation
    'il_4_concentration_picogram_per_milliliter',    # Th2, B-cell activation
    'il_10_concentration_picogram_per_milliliter',   # anti-inflammatory, regulatory
    'il_13_concentration_picogram_per_milliliter',   # Th2 effector
    'il_15_concentration_picogram_per_milliliter',   # NK/T-cell homeostasis
    'il_17a_concentration_picogram_per_milliliter',  # Th17 axis
    'il_1ra_concentration_picogram_per_milliliter',  # IL-1 antagonist, regulatory brake
    'il_21_concentration_picogram_per_milliliter',   # T follicular helper
    'il_27_concentration_picogram_per_milliliter',   # Treg induction
    'tnfb_concentration_picogram_per_milliliter',    # immunosuppressive (if this is TGF-α)
    'ifny_concentration_picogram_per_milliliter',    # Th1 effector
]

reg_df = immune_eve_long[
    immune_eve_long['timepoint'].str.startswith('R+') &
    immune_eve_long['marker'].isin(immune_regulation_markers)
]

composite = (
    reg_df.groupby(['astronaut', 'timepoint'])['marker_score'] 
    .mean()
    .reset_index()
    .rename(columns={'marker_score': 'immune_regulation_score'})
)
composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
composite = composite.sort_values(['astronaut', 'timepoint_num'])

overall_scores = (
    composite
    .groupby('astronaut')['immune_regulation_score']
    .mean()
    .reset_index()
    .rename(columns={'immune_regulation_score': 'overall_immune_regulation_score'})
    .sort_values('overall_immune_regulation_score', ascending=False)
)

print(overall_scores)
