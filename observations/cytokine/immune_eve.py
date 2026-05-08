from cytokine_dataframe import *

# =============================================================================
# LOAD & PROCESS
# =============================================================================

immune_eve = organize_df(
    pd.read_csv(
        'https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download'
        '?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv',
        index_col=0
    ).transpose()
)

immune_eve_long = process_cytokine_df(immune_eve, 'immune_eve')

# =============================================================================
# DOMAIN MARKER LISTS
# Full marker names exactly as they appear in the dataset.
# Comments indicate the biological role within the domain.
# =============================================================================

MARKERS = {
    'immune_regulation': [
        'il_2_concentration_picogram_per_milliliter',    # T-cell proliferation
        'il_4_concentration_picogram_per_milliliter',    # Th2, B-cell activation
        'il_10_concentration_picogram_per_milliliter',   # anti-inflammatory, regulatory brake
        'il_13_concentration_picogram_per_milliliter',   # Th2 effector
        'il_15_concentration_picogram_per_milliliter',   # NK/T-cell homeostasis
        'il_17a_concentration_picogram_per_milliliter',  # Th17 axis
        'il_1ra_concentration_picogram_per_milliliter',  # IL-1 antagonist, regulatory brake
        'il_21_concentration_picogram_per_milliliter',   # T follicular helper
        'il_27_concentration_picogram_per_milliliter',   # Treg induction
        'ifny_concentration_picogram_per_milliliter',    # Th1 effector
    ],
    'inflammation': [
        'il_1a_concentration_picogram_per_milliliter',   # acute phase inflammatory
        'il_1b_concentration_picogram_per_milliliter',   # canonical inflammasome cytokine
        'il_6_concentration_picogram_per_milliliter',    # acute phase, fever response
        'il_8_concentration_picogram_per_milliliter',    # neutrophil chemoattractant
        'il_18_concentration_picogram_per_milliliter',   # inflammasome, IFN-γ inducer
        'il_33_concentration_picogram_per_milliliter',   # alarmin, tissue inflammation
        'tnfa_concentration_picogram_per_milliliter',    # master pro-inflammatory cytokine
        'tnfb_concentration_picogram_per_milliliter',    # pro-inflammatory
        'groa_concentration_picogram_per_milliliter',    # neutrophil chemoattractant (CXCL1)
        'mcp_1_concentration_picogram_per_milliliter',   # monocyte recruitment
        'mip_1a_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
        'mip_1b_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
        'ip_10_concentration_picogram_per_milliliter',   # CXCL10, IFN-driven inflammation
    ],
}

# =============================================================================
# MARKER WEIGHTS
# 1.0 = primary driver of the domain (well-established central role)
# 0.5 = secondary contributor (associated but not a master regulator)
# =============================================================================

MARKER_WEIGHTS = {
    # Immune regulation
    'il_2_concentration_picogram_per_milliliter':   1.0,  # central T-cell signal
    'il_4_concentration_picogram_per_milliliter':   0.5,  # secondary Th2
    'il_10_concentration_picogram_per_milliliter':  1.0,  # key regulatory cytokine
    'il_13_concentration_picogram_per_milliliter':  0.5,  # secondary Th2 effector
    'il_15_concentration_picogram_per_milliliter':  1.0,  # NK/T homeostasis
    'il_17a_concentration_picogram_per_milliliter': 0.5,  # secondary Th17
    'il_1ra_concentration_picogram_per_milliliter': 1.0,  # important IL-1 regulator
    'il_21_concentration_picogram_per_milliliter':  0.5,  # secondary T follicular
    'il_27_concentration_picogram_per_milliliter':  0.5,  # secondary Treg
    'ifny_concentration_picogram_per_milliliter':   1.0,  # Th1 master regulator

    # Inflammation
    'il_1a_concentration_picogram_per_milliliter':  0.5,  # secondary acute phase
    'il_1b_concentration_picogram_per_milliliter':  1.0,  # canonical inflammasome
    'il_6_concentration_picogram_per_milliliter':   1.0,  # master acute phase cytokine
    'il_8_concentration_picogram_per_milliliter':   1.0,  # primary neutrophil signal
    'il_18_concentration_picogram_per_milliliter':  0.5,  # secondary inflammasome
    'il_33_concentration_picogram_per_milliliter':  0.5,  # secondary alarmin
    'tnfa_concentration_picogram_per_milliliter':   1.0,  # master pro-inflammatory
    'tnfb_concentration_picogram_per_milliliter':   0.5,  # secondary pro-inflammatory
    'groa_concentration_picogram_per_milliliter':   0.5,  # secondary chemoattractant
    'mcp_1_concentration_picogram_per_milliliter':  1.0,  # primary monocyte recruiter
    'mip_1a_concentration_picogram_per_milliliter': 0.5,  # secondary MIP
    'mip_1b_concentration_picogram_per_milliliter': 0.5,  # secondary MIP
    'ip_10_concentration_picogram_per_milliliter':  1.0,  # primary IFN-driven marker
}

# Register weights before scoring
set_marker_weights(MARKER_WEIGHTS)

# =============================================================================
# SCORE
# =============================================================================

print("Immune Regulation Score:")
immune_reg_tp, immune_reg_overall = score(
    immune_eve_long,
    MARKERS['immune_regulation'],
    domain_name='Immune Regulation'
)

print("\n\nInflammation Score:")
inflammation_tp, inflammation_overall = score(
    immune_eve_long,
    MARKERS['inflammation'],
    domain_name='Inflammation'
)

# from cytokine_dataframe import *

# # Dataset straight from open-source
# immune_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv', index_col=0).transpose())

# # Clean up dataset
# immune_eve_long = process_cytokine_df(immune_eve, 'immune_eve')

# # # Make observations on significant changes in concentration before and after flight
# # postflight = immune_eve_long[immune_eve_long['timepoint'].str.startswith('R+')]
# # significant = postflight[postflight['marker_score'].abs() >= 2]

# # # Make a dataframe of the persistence of significant changes in certain markers over timepoints
# # persistence = significant.groupby(['astronaut', 'marker']).agg(
# #     times_significant = ('marker_score', 'count'), # count rows
# #     max_score = ('marker_score', 'max'), # highest score
# #     mean_log2fc = ('log2fc', 'mean') # average log2fc
# # ).reset_index()

# ### BEGIN CLASSIFICATION ###

# markers = {
#     'immune_regulation': [
#     'il_2_concentration_picogram_per_milliliter',    # T-cell proliferation
#     'il_4_concentration_picogram_per_milliliter',    # Th2, B-cell activation
#     'il_10_concentration_picogram_per_milliliter',   # anti-inflammatory, regulatory
#     'il_13_concentration_picogram_per_milliliter',   # Th2 effector
#     'il_15_concentration_picogram_per_milliliter',   # NK/T-cell homeostasis
#     'il_17a_concentration_picogram_per_milliliter',  # Th17 axis
#     'il_1ra_concentration_picogram_per_milliliter',  # IL-1 antagonist, regulatory brake
#     'il_21_concentration_picogram_per_milliliter',   # T follicular helper
#     'il_27_concentration_picogram_per_milliliter',   # Treg induction
#     'ifny_concentration_picogram_per_milliliter'     # Th1 effector
#     ],'inflammation': [
#     'il_1a_concentration_picogram_per_milliliter',   # acute phase inflammatory
#     'il_1b_concentration_picogram_per_milliliter',   # canonical inflammasome cytokine
#     'il_6_concentration_picogram_per_milliliter',    # acute phase, fever response
#     'il_8_concentration_picogram_per_milliliter',    # neutrophil chemoattractant
#     'il_18_concentration_picogram_per_milliliter',   # inflammasome, IFN-γ inducer
#     'il_33_concentration_picogram_per_milliliter',   # alarmin, tissue inflammation
#     'tnfa_concentration_picogram_per_milliliter',    # master pro-inflammatory cytokine
#     'tnfb_concentration_picogram_per_milliliter',    # pro-inflammatory
#     'groa_concentration_picogram_per_milliliter',    # neutrophil chemoattractant (CXCL1)
#     'mcp_1_concentration_picogram_per_milliliter',   # monocyte recruitment
#     'mip_1a_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
#     'mip_1b_concentration_picogram_per_milliliter',  # macrophage inflammatory protein
#     'ip_10_concentration_picogram_per_milliliter'    # CXCL10, IFN-driven inflammation
# ]}

# # print("Immune Regulation Score:")
# # immune_reg_composite, immune_reg_overall_score = score(immune_eve_long, markers['immune_regulation'])
# # print("\nInflammation Score:")
# # inflammation_composite, imflammation_overall_score = score(immune_eve_long, markers['inflammation'])
# # to_plot = [immune_reg_composite, inflammation_composite]
# immune_reg_tp, immune_reg_overall = score(immune_eve_long, markers['immune_regulation'], 'Immune Regulation')
# inflammation_tp, inflammation_overall = score(immune_eve_long, markers['inflammation'], 'Inflammation')
