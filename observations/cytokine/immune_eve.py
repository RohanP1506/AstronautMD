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

markers = {
    'immune_regulation': [
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
    'ifny_concentration_picogram_per_milliliter'     # Th1 effector
    ],'inflammation': [
    'il_1a_concentration_picogram_per_milliliter',   # acute phase inflammatory
    'il_1b_concentration_picogram_per_milliliter',   # canonical inflammasome cytokine
    'il_6_concentration_picogram_per_milliliter',    # acute phase, fever response
    'il_8_concentration_picogram_per_milliliter',    # neutrophil chemoattractant
    'il_18_concentration_picogram_per_milliliter',   # inflammasome, IFN-γ inducer
    'il_33_concentration_picogram_per_milliliter',   # alarmin, tissue inflammation
    'tnfa_concentration_picogram_per_milliliter',    # master pro-inflammatory cytokine
    'ifny_concentration_picogram_per_milliliter',    # Th1 inflammation (shared with immune reg)
    'groa_concentration_picogram_per_milliliter',    # neutrophil chemoattractant (CXCL1)
    'mcp_1_concentration_picogram_per_milliliter',   # monocyte recruitment
    'mip_1a_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
    'mip_1b_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
    'ip_10_concentration_picogram_per_milliliter'    # CXCL10, IFN-driven inflammation
]}

print("Immune Regulation Score:")
immune_reg_composite, immune_reg_overall_score = score(immune_eve_long, markers['immune_regulation'])
print("\nInflammation Score:")
inflammation_composite, imflammation_overall_score = score(immune_eve_long, markers['inflammation'])
to_plot = [immune_reg_composite, inflammation_composite]
