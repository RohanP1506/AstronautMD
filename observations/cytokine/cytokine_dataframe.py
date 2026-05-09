import numpy as np
import pandas as pd

# =============================================================================
# DATA LOADING & PROCESSING
# =============================================================================

def organize_df(df):
    """
    Reshape a transposed NASA multiplex serum CSV into long format.
    Only keeps concentration rows. Drops NaN and zero concentrations.
    Output columns: marker, sample_id, concentration, astronaut, timepoint
    """
    conc = df.reset_index()
    conc = conc[conc['index'].str.contains('concentration')]
    conc = conc.rename(columns={'index': 'marker'})
    melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')
    melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')
    melted_df = melted_df.dropna(subset=['astronaut', 'timepoint'])
    melted_df['concentration'] = pd.to_numeric(melted_df['concentration'], errors='coerce')
    melted_df = melted_df[melted_df['concentration'] > 0].dropna(subset=['concentration'])
    return melted_df.reset_index(drop=True)


def compute_baseline(df):
    """
    Personal baseline per astronaut per marker: mean across all L- preflight timepoints.
    """
    preflight = df[df['timepoint'].str.startswith('L-')]
    return (
        preflight
        .groupby(['astronaut', 'marker'])['concentration']
        .mean()
        .reset_index()
        .rename(columns={'concentration': 'baseline_mean'})
    )


def process_cytokine_df(df, name):
    """
    Merge baseline and compute log2FC. No scoring here.
    Guards against log2(0) or negative concentrations.
    """
    baseline = compute_baseline(df)
    df = df.merge(baseline, on=['astronaut', 'marker'])
    df = df[(df['concentration'] > 0) & (df['baseline_mean'] > 0)]
    df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    df['source'] = name
    return df.reset_index(drop=True)


# =============================================================================
# SCORING CONFIGURATION
#
# Methodology:
#   - Significance threshold: |log2FC| >= 1.0  (2x fold-change)
#   - Domain score at each timepoint = weighted fraction of significant markers
#       weighted_fraction = sum(weights of significant markers) /
#                           sum(weights of all domain markers)
#   - Marker weights: 1.0 = primary contributor, 0.5 = secondary contributor
#   - Overall domain score = timepoint-weighted mean of timepoint fractions
#   - Timepoint weights: R+1 down-weighted to 0.5 (reentry physical stress confound)
#   - Risk tier: overall fraction binned to 1-5 scale
# =============================================================================

SIGNIFICANCE_THRESHOLD = 1.0  # |log2FC| >= 1.0 = 2x fold-change

TIMEPOINT_WEIGHTS = {
    'R+1':   0.5,   # down-weighted — physical reentry stress confound
    'R+45':  1.0,
    'R+82':  1.0,
    'R+194': 1.0,
}

_MARKER_WEIGHTS = {}

def set_marker_weights(weights):
    """Register marker importance weights before calling score()."""
    global _MARKER_WEIGHTS
    _MARKER_WEIGHTS = weights


def _fraction_to_risk(fraction):
    """Bin weighted fraction of significant markers into 1-5 risk score."""
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                 return 5


def _clean_name(marker):
    return marker.split('_concentration')[0].upper()


# =============================================================================
# MAIN SCORING FUNCTION
# =============================================================================

def score(df, markers, domain_name=''):
    """
    Score a biological domain for each astronaut using postflight serum data.

    Returns
    -------
    timepoint_df : per astronaut x timepoint summary
    overall_df   : per astronaut overall weighted fraction and 1-5 risk score
    """
    domain_df = df[
        df['timepoint'].str.startswith('R+') &
        df['marker'].isin(markers)
    ].copy()

    if domain_df.empty:
        print(f"Warning: No data found for domain '{domain_name}'. Check marker names.")
        return pd.DataFrame(), pd.DataFrame()

    domain_df['significant']   = domain_df['log2fc'].abs() >= SIGNIFICANCE_THRESHOLD
    domain_df['direction']     = np.where(domain_df['log2fc'] > 0, 'up', 'down')
    domain_df['marker_weight'] = domain_df['marker'].map(_MARKER_WEIGHTS).fillna(0.5)

    def summarize_timepoint(group):
        sig          = group[group['significant']]
        total_weight = group['marker_weight'].sum()
        weighted_sig = sig['marker_weight'].sum()
        fraction     = weighted_sig / total_weight if total_weight > 0 else 0.0

        up_df   = sig[sig['direction'] == 'up'][['marker', 'log2fc']].copy()
        down_df = sig[sig['direction'] == 'down'][['marker', 'log2fc']].copy()
        up_df['marker']   = up_df['marker'].apply(_clean_name)
        down_df['marker'] = down_df['marker'].apply(_clean_name)

        return pd.Series({
            'fraction':    round(fraction, 4),
            'n_sig':       len(sig),
            'n_total':     len(group),
            'up_scores':   dict(zip(up_df['marker'],   up_df['log2fc'].round(2))),
            'down_scores': dict(zip(down_df['marker'], down_df['log2fc'].round(2))),
        })

    timepoint_df = (
        domain_df
        .groupby(['astronaut', 'timepoint'])
        .apply(summarize_timepoint)
        .reset_index()
    )

    timepoint_df['_tp_num']               = timepoint_df['timepoint'].str.extract(r'R\+(\d+)').astype(int)
    timepoint_df['timepoint_weight']      = timepoint_df['timepoint'].map(TIMEPOINT_WEIGHTS).fillna(1.0)
    timepoint_df['weighted_contribution'] = timepoint_df['fraction'] * timepoint_df['timepoint_weight']
    timepoint_df = timepoint_df.sort_values(['astronaut', '_tp_num']).drop(columns='_tp_num')

    def overall_score(g):
        return pd.Series({
            'overall_fraction': g['weighted_contribution'].sum() / g['timepoint_weight'].sum()
        })

    overall_df = (
        timepoint_df
        .groupby('astronaut')
        .apply(overall_score)
        .reset_index()
    )
    overall_df['risk_score'] = overall_df['overall_fraction'].apply(_fraction_to_risk)
    overall_df = overall_df.sort_values('risk_score', ascending=False).reset_index(drop=True)

    # Printout
    print(f"\n{'='*65}")
    print(f"  {domain_name}")
    print(f"  Threshold: |log2FC| >= {SIGNIFICANCE_THRESHOLD} (2x fold-change)")
    print(f"  R+1 weighted 0.5x (reentry stress confound)")
    print(f"{'='*65}")

    for astronaut, ast_df in timepoint_df.groupby('astronaut'):
        risk = overall_df.loc[overall_df['astronaut'] == astronaut, 'risk_score'].values[0]
        frac = overall_df.loc[overall_df['astronaut'] == astronaut, 'overall_fraction'].values[0]
        print(f"\n  {astronaut}  —  Risk score: {risk}/5  (weighted fraction: {frac:.2f})")

        for _, row in ast_df.iterrows():
            print(f"\n    {row['timepoint']}  ({int(row['n_sig'])}/{int(row['n_total'])} markers significant)")

            if row['up_scores']:
                print(f"      Elevated:")
                for marker, fc in sorted(row['up_scores'].items(), key=lambda x: -x[1]):
                    print(f"        + {marker:<25} log2FC: {fc:+.2f}  ({2**fc:.1f}x increase)")

            if row['down_scores']:
                print(f"      Suppressed:")
                for marker, fc in sorted(row['down_scores'].items(), key=lambda x: x[1]):
                    print(f"        - {marker:<25} log2FC: {fc:+.2f}  ({2**abs(fc):.1f}x decrease)")

            if not row['up_scores'] and not row['down_scores']:
                print(f"        No significant markers at this timepoint")

    print(f"\n{'─'*65}")
    print(f"  Overall risk scores:")
    print(overall_df[['astronaut', 'overall_fraction', 'risk_score']].round(3).to_string(index=False))

    return timepoint_df, overall_df

# import numpy as np
# import pandas as pd

# # =============================================================================
# # DATA LOADING & PROCESSING
# # =============================================================================

# def organize_df(df):
#     """
#     Reshape a transposed NASA multiplex CSV into long format.
#     Input:  transposed CSV where rows = samples, columns = markers
#     Output: long dataframe with columns: marker, sample_id, concentration,
#             astronaut, timepoint
#     Only keeps concentration columns (drops percent columns).
#     """
#     conc = df.reset_index()
#     conc = conc[conc['index'].str.contains('concentration')]
#     conc = conc.rename(columns={'index': 'marker'})
#     melted_df = conc.melt(id_vars='marker', var_name='sample_id', value_name='concentration')
#     melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')
#     melted_df = melted_df.dropna(subset=['astronaut', 'timepoint'])
#     melted_df['concentration'] = pd.to_numeric(melted_df['concentration'], errors='coerce')
#     melted_df = melted_df.dropna(subset=['concentration'])
#     return melted_df.reset_index(drop=True)


# def compute_baseline(df):
#     """
#     Compute each astronaut's personal baseline per marker as the mean
#     concentration across all preflight timepoints (L-92, L-44, L-3).
#     """
#     preflight = df[df['timepoint'].str.startswith('L-')]
#     return (
#         preflight
#         .groupby(['astronaut', 'marker'])['concentration']
#         .mean()
#         .reset_index()
#         .rename(columns={'concentration': 'baseline_mean'})
#     )


# def process_cytokine_df(df, name):
#     """
#     Merge baseline into the dataframe and compute log2FC per marker per
#     timepoint relative to that astronaut's personal preflight baseline.
#     No scoring is done here — scoring happens in score().
#     """
#     baseline = compute_baseline(df)
#     df = df.merge(baseline, on=['astronaut', 'marker'])

#     # Guard against log2(0) or negative concentrations
#     df = df[(df['concentration'] > 0) & (df['baseline_mean'] > 0)]

#     df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
#     df['source'] = name
#     return df.reset_index(drop=True)


# # =============================================================================
# # SCORING CONFIGURATION
# #
# # Methodology:
# #   - Significance threshold: |log2FC| >= 1.0  (2x fold-change)
# #   - Domain score at each timepoint = weighted fraction of significant markers
# #       weighted_fraction = sum(weights of significant markers) /
# #                           sum(weights of all domain markers)
# #   - Marker weights: 1.0 = primary contributor, 0.5 = secondary contributor
# #   - Overall domain score = timepoint-weighted mean of timepoint fractions
# #   - Timepoint weights: R+1 down-weighted to 0.5 (reentry stress confound);
# #       R+45, R+82, R+194 weighted 1.0
# #   - Risk tier: overall fraction binned to 1-5 scale
# # =============================================================================

# SIGNIFICANCE_THRESHOLD = 1.0   # |log2FC| >= 1.0 = 2x fold-change

# TIMEPOINT_WEIGHTS = {
#     'R+1':   0.5,   # down-weighted — physical reentry stress is a confound
#     'R+45':  1.0,
#     'R+82':  1.0,
#     'R+194': 1.0,
# }

# # Marker weights are set per-panel by calling set_marker_weights()
# _MARKER_WEIGHTS = {}

# def set_marker_weights(weights):
#     """
#     Register marker importance weights before calling score().
#     weights: dict mapping full marker name -> float (1.0 or 0.5)
#     Any marker not in this dict defaults to 0.5.
#     """
#     global _MARKER_WEIGHTS
#     _MARKER_WEIGHTS = weights


# def _fraction_to_risk(fraction):
#     """
#     Bin a weighted fraction of significant markers into a 1-5 risk score.
#       1 = minimal    < 10% of weighted markers affected
#       2 = mild       10-25% affected
#       3 = moderate   25-50% affected
#       4 = significant 50-75% affected
#       5 = severe     > 75% affected
#     """
#     if fraction < 0.10:    return 1
#     elif fraction < 0.25:  return 2
#     elif fraction < 0.50:  return 3
#     elif fraction < 0.75:  return 4
#     else:                  return 5


# def _clean_name(marker):
#     """
#     Strip unit suffix for readable output.
#     e.g. 'il_6_concentration_picogram_per_milliliter' -> 'IL_6'
#     """
#     return marker.split('_concentration')[0].upper()


# # =============================================================================
# # MAIN SCORING FUNCTION
# # =============================================================================

# def score(df, markers, domain_name=''):
#     """
#     Score a biological domain for each astronaut using postflight data.

#     Parameters
#     ----------
#     df          : processed dataframe from process_cytokine_df() with log2fc column
#     markers     : list of full marker names belonging to this domain
#     domain_name : label used in printouts

#     Returns
#     -------
#     timepoint_df : per astronaut x timepoint summary with fraction, n_sig,
#                    n_total, up_scores, down_scores
#     overall_df   : per astronaut overall weighted fraction and 1-5 risk score
#     """
#     # Filter to postflight timepoints and domain markers only
#     domain_df = df[
#         df['timepoint'].str.startswith('R+') &
#         df['marker'].isin(markers)
#     ].copy()

#     if domain_df.empty:
#         print(f"⚠️  No data found for domain '{domain_name}'. Check marker names.")
#         return pd.DataFrame(), pd.DataFrame()

#     # Flag significance and direction
#     domain_df['significant']    = domain_df['log2fc'].abs() >= SIGNIFICANCE_THRESHOLD
#     domain_df['direction']      = np.where(domain_df['log2fc'] > 0, 'up', 'down')
#     domain_df['marker_weight']  = domain_df['marker'].map(_MARKER_WEIGHTS).fillna(0.5)

#     # ── Per-timepoint summary ──────────────────────────────────────────────────
#     def summarize_timepoint(group):
#         sig          = group[group['significant']]
#         total_weight = group['marker_weight'].sum()
#         weighted_sig = sig['marker_weight'].sum()
#         fraction     = weighted_sig / total_weight if total_weight > 0 else 0.0

#         up_df   = sig[sig['direction'] == 'up'][['marker', 'log2fc']].copy()
#         down_df = sig[sig['direction'] == 'down'][['marker', 'log2fc']].copy()
#         up_df['marker']   = up_df['marker'].apply(_clean_name)
#         down_df['marker'] = down_df['marker'].apply(_clean_name)

#         return pd.Series({
#             'fraction':    round(fraction, 4),
#             'n_sig':       len(sig),
#             'n_total':     len(group),
#             'up_scores':   dict(zip(up_df['marker'],   up_df['log2fc'].round(2))),
#             'down_scores': dict(zip(down_df['marker'], down_df['log2fc'].round(2))),
#         })

#     timepoint_df = (
#         domain_df
#         .groupby(['astronaut', 'timepoint'])
#         .apply(summarize_timepoint)
#         .reset_index()
#     )

#     # Sort timepoints numerically (R+1, R+45, R+82, R+194)
#     timepoint_df['_tp_num']             = timepoint_df['timepoint'].str.extract(r'R\+(\d+)').astype(int)
#     timepoint_df['timepoint_weight']    = timepoint_df['timepoint'].map(TIMEPOINT_WEIGHTS).fillna(1.0)
#     timepoint_df['weighted_contribution'] = timepoint_df['fraction'] * timepoint_df['timepoint_weight']
#     timepoint_df = timepoint_df.sort_values(['astronaut', '_tp_num']).drop(columns='_tp_num')

#     # ── Overall score per astronaut ────────────────────────────────────────────
#     def overall_score(g):
#         return pd.Series({
#             'overall_fraction': g['weighted_contribution'].sum() / g['timepoint_weight'].sum()
#         })

#     overall_df = (
#         timepoint_df
#         .groupby('astronaut')
#         .apply(overall_score)
#         .reset_index()
#     )
#     overall_df['risk_score'] = overall_df['overall_fraction'].apply(_fraction_to_risk)
#     overall_df = overall_df.sort_values('risk_score', ascending=False).reset_index(drop=True)

#     # ── Printout ───────────────────────────────────────────────────────────────
#     print(f"\n{'='*65}")
#     print(f"  {domain_name}")
#     print(f"  Significance: |log2FC| >= {SIGNIFICANCE_THRESHOLD} (2x fold-change)")
#     print(f"  R+1 weighted 0.5x (reentry stress confound)")
#     print(f"{'='*65}")

#     for astronaut, ast_df in timepoint_df.groupby('astronaut'):
#         risk = overall_df.loc[overall_df['astronaut'] == astronaut, 'risk_score'].values[0]
#         frac = overall_df.loc[overall_df['astronaut'] == astronaut, 'overall_fraction'].values[0]
#         print(f"\n  {astronaut}  —  Risk score: {risk}/5  (weighted fraction: {frac:.2f})")

#         for _, row in ast_df.iterrows():
#             print(f"\n    {row['timepoint']}  "
#                   f"({int(row['n_sig'])}/{int(row['n_total'])} markers significant)")

#             if row['up_scores']:
#                 print(f"      Elevated:")
#                 for marker, fc in sorted(row['up_scores'].items(), key=lambda x: -x[1]):
#                     print(f"        + {marker:<25} log2FC: {fc:+.2f}  ({2**fc:.1f}x increase)")

#             if row['down_scores']:
#                 print(f"      Suppressed:")
#                 for marker, fc in sorted(row['down_scores'].items(), key=lambda x: x[1]):
#                     print(f"        - {marker:<25} log2FC: {fc:+.2f}  ({2**abs(fc):.1f}x decrease)")

#             if not row['up_scores'] and not row['down_scores']:
#                 print(f"        No significant markers at this timepoint")

#     print(f"\n{'─'*65}")
#     print(f"  Overall risk scores:")
#     print(overall_df[['astronaut', 'overall_fraction', 'risk_score']].round(3).to_string(index=False))

#     return timepoint_df, overall_df


# # import numpy as np
# # import pandas as pd

# # # Clean dataframe for easy studies
# # def organize_df(df):
# #     # Modify dataframe to focus only on concentrations
# #     conc = df.reset_index()
# #     conc = conc[conc['index'].str.contains('concentration')]
# #     conc = conc.rename(columns={'index': 'marker'})

# #     # Melt dataframe into 5 different columns: marker, sample id, concentration, astronaut and timepoint
# #     melted_df = conc.melt(id_vars = 'marker', var_name = 'sample_id', value_name = 'concentration')
# #     melted_df[['astronaut', 'timepoint']] = melted_df['sample_id'].str.extract(r'(C\d+)_serum_(.*)')

# #     return melted_df

# # # Get the baseline mean
# # def compute_baseline(df):
# #     # Find whatever preflight timepoints exist in this dataset
# #     preflight_mask = df['timepoint'].str.startswith('L-')
# #     preflight = df[preflight_mask]
    
# #     # Average them per astronaut per marker
# #     baseline = preflight.groupby(['astronaut', 'marker'])['concentration'].mean().reset_index()
# #     baseline = baseline.rename(columns={'concentration': 'baseline_mean'})

# #     return baseline












# # def process_cytokine_df(df, name):
# #     # Compute personal baseline from whatever preflight timepoints exist
# #     baseline = compute_baseline(df)
    
# #     # Merge and compute log2fc
# #     df = df.merge(baseline, on=['astronaut', 'marker'])
# #     df['log2fc'] = np.log2(df['concentration'] / df['baseline_mean'])
    
# #     # Score
# #     df['marker_score'] = df['log2fc'].apply(score_log2fc)
    
# #     # Tag which dataset this came from
# #     df['source'] = name
    
# #     return df

# # def score(df, markers):
# #     df = df[
# #         df['timepoint'].str.startswith('R+') &
# #         df['marker'].isin(markers)
# #     ]

# #     composite = (
# #         df.groupby(['astronaut', 'timepoint'])['marker_score']
# #         .mean()
# #         .reset_index()
# #     )
# #     composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
# #     composite = composite.sort_values(['astronaut', 'timepoint_num'])

# #     # Print per-timepoint scores per astronaut
# #     print("\nPer-timepoint scores:")
# #     per_tp = (
# #         composite
# #         .pivot_table(index='astronaut', columns='timepoint', values='marker_score')
# #     )
# #     ordered_cols = sorted(per_tp.columns, key=lambda x: int(x.replace('R+','')))
# #     per_tp = per_tp[ordered_cols].round(2)
# #     print(per_tp.to_string())

# #     overall_scores = (
# #         composite
# #         .groupby('astronaut')['marker_score']
# #         .mean()
# #         .reset_index()
# #         .rename(columns={'marker_score': 'overall_score'})
# #         .sort_values('overall_score', ascending=False)
# #     )

# #     # Print overall weighted scores
# #     print("\nOverall weighted score:")
# #     print(overall_scores.round(2).to_string(index=False))

# #     return composite, overall_scores
