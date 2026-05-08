import numpy as np
import pandas as pd

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
#   - Timepoint weights: R+1 down-weighted to 0.5 (reentry stress confound)
#   - Risk tier: overall fraction binned to 1-5 scale
# =============================================================================

SIGNIFICANCE_THRESHOLD = 1.0   # |log2FC| >= 1.0 = 2x fold-change

TIMEPOINT_WEIGHTS = {
    'R+1':   0.5,   # down-weighted — physical reentry stress is a confound
    'R+45':  1.0,
    'R+82':  1.0,
    'R+194': 1.0,
}

# Domain marker lists for the urine OSD-656 Alamar Panel
# Keys are short substrings that will match inside full marker names
DOMAIN_MAP = {
    'inflammation': [
        'il6', 'il1b', 'tnf', 'ccl2', 'ccl25', 'cxcl2', 'mmp3', 'mmp8', 'ptx3'
    ],
    'oxidative_stress': [
        'ager', 'gzmb', '8ohdg', 'hgf'
    ]
}

# Marker importance weights — matched by substring against full marker names
# 1.0 = primary driver, 0.5 = secondary contributor
DOMAIN_WEIGHTS = {
    # Inflammation
    'il6':   1.0,   # master acute phase cytokine
    'il1b':  1.0,   # canonical inflammasome cytokine
    'tnf':   1.0,   # master pro-inflammatory
    'ccl2':  1.0,   # primary monocyte recruiter
    'ccl25': 0.5,   # secondary chemokine
    'cxcl2': 0.5,   # secondary neutrophil chemoattractant
    'mmp3':  0.5,   # secondary tissue remodeling
    'mmp8':  0.5,   # secondary tissue remodeling
    'ptx3':  1.0,   # acute phase protein, primary inflammation marker

    # Oxidative stress
    'ager':  1.0,   # RAGE — primary oxidative stress receptor
    'gzmb':  0.5,   # secondary cytotoxic marker
    '8ohdg': 1.0,   # 8-OHdG — primary DNA oxidation marker
    'hgf':   0.5,   # secondary growth/repair factor
}


# =============================================================================
# HELPERS
# =============================================================================

def _find_domain(marker_name, domain_mapping):
    """Return the domain for a marker by substring match, or None."""
    marker_clean = str(marker_name).lower()
    for domain, keywords in domain_mapping.items():
        for kw in keywords:
            if kw.lower() in marker_clean:
                return domain
    return None


def _get_marker_weight(marker_name, weight_map):
    """Return the importance weight for a marker by substring match."""
    marker_clean = str(marker_name).lower()
    for kw, weight in weight_map.items():
        if kw.lower() in marker_clean:
            return weight
    return 0.5   # default to secondary if not found


def _clean_name(marker):
    """Strip unit suffix for readable output."""
    return marker.split('_concentration')[0].upper()


def _fraction_to_risk(fraction):
    """
    Bin a weighted fraction of significant markers into a 1-5 risk score.
      1 = minimal     < 10% of weighted markers affected
      2 = mild        10-25% affected
      3 = moderate    25-50% affected
      4 = significant 50-75% affected
      5 = severe      > 75% affected
    """
    if fraction < 0.10:    return 1
    elif fraction < 0.25:  return 2
    elif fraction < 0.50:  return 3
    elif fraction < 0.75:  return 4
    else:                  return 5


# =============================================================================
# MAIN SCORING FUNCTION
# =============================================================================

def calculate_composite_scores(processed_df, domain_mapping, weight_map=None, domain_name=''):
    """
    Score biological domains for each astronaut using postflight urine data.

    Parameters
    ----------
    processed_df  : dataframe from process_urine_df() with log2fc column
    domain_mapping: dict of {domain: [keyword substrings]}
    weight_map    : dict of {keyword: float} for marker importance weights.
                    Defaults to DOMAIN_WEIGHTS if not provided.
    domain_name   : optional label for printouts

    Returns
    -------
    timepoint_df : per astronaut x timepoint summary
    overall_df   : per astronaut overall weighted fraction and 1-5 risk score
    """
    if weight_map is None:
        weight_map = DOMAIN_WEIGHTS

    # Assign domain and weight to each marker row
    df = processed_df.copy()
    df['domain']         = df['marker'].apply(lambda m: _find_domain(m, domain_mapping))
    df['marker_weight']  = df['marker'].apply(lambda m: _get_marker_weight(m, weight_map))

    # Drop markers not assigned to any domain
    df = df.dropna(subset=['domain'])

    if df.empty:
        print("⚠️  No domain matches found. Check that marker names contain the keywords in DOMAIN_MAP.")
        print("    Sample marker names:", processed_df['marker'].unique()[:5].tolist())
        return pd.DataFrame(), pd.DataFrame()

    # Filter to postflight timepoints only
    df = df[df['timepoint'].str.startswith('R+')].copy()

    # Flag significance and direction
    df['significant'] = df['log2fc'].abs() >= SIGNIFICANCE_THRESHOLD
    df['direction']   = np.where(df['log2fc'] > 0, 'up', 'down')

    # ── Per-timepoint summary ──────────────────────────────────────────────────
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
        df
        .groupby(['astronaut', 'timepoint', 'domain'])
        .apply(summarize_timepoint)
        .reset_index()
    )

    # Sort timepoints numerically
    timepoint_df['_tp_num']               = timepoint_df['timepoint'].str.extract(r'R\+(\d+)').astype(int)
    timepoint_df['timepoint_weight']      = timepoint_df['timepoint'].map(TIMEPOINT_WEIGHTS).fillna(1.0)
    timepoint_df['weighted_contribution'] = timepoint_df['fraction'] * timepoint_df['timepoint_weight']
    timepoint_df = timepoint_df.sort_values(['astronaut', 'domain', '_tp_num']).drop(columns='_tp_num')

    # ── Overall score per astronaut per domain ─────────────────────────────────
    overall_df = (
        timepoint_df
        .groupby(['astronaut', 'domain'])
        .apply(lambda g: pd.Series({
            'overall_fraction': g['weighted_contribution'].sum() / g['timepoint_weight'].sum()
        }))
        .reset_index()
    )
    overall_df['risk_score'] = overall_df['overall_fraction'].apply(_fraction_to_risk)
    overall_df = overall_df.sort_values(['domain', 'risk_score'], ascending=[True, False]).reset_index(drop=True)

    # ── Printout ───────────────────────────────────────────────────────────────
    for domain in timepoint_df['domain'].unique():
        print(f"\n{'='*65}")
        print(f"  {domain.replace('_', ' ').title()}")
        print(f"  Significance: |log2FC| >= {SIGNIFICANCE_THRESHOLD} (2x fold-change)")
        print(f"  R+1 weighted 0.5x (reentry stress confound)")
        print(f"{'='*65}")

        domain_tp  = timepoint_df[timepoint_df['domain'] == domain]
        domain_ov  = overall_df[overall_df['domain'] == domain]

        for astronaut, ast_df in domain_tp.groupby('astronaut'):
            risk = domain_ov.loc[domain_ov['astronaut'] == astronaut, 'risk_score'].values[0]
            frac = domain_ov.loc[domain_ov['astronaut'] == astronaut, 'overall_fraction'].values[0]
            print(f"\n  {astronaut}  —  Risk score: {risk}/5  (weighted fraction: {frac:.2f})")

            for _, row in ast_df.iterrows():
                print(f"\n    {row['timepoint']}  "
                      f"({int(row['n_sig'])}/{int(row['n_total'])} markers significant)")

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
        print(f"  Overall risk scores ({domain}):")
        print(domain_ov[['astronaut', 'overall_fraction', 'risk_score']].round(3).to_string(index=False))

    return timepoint_df, overall_df