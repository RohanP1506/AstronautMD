"""
cbc_score.py — CBC domain scoring for AstronautMD
Scores CBC log2FC data per astronaut per canonical biological domain.
"""

import numpy as np
import pandas as pd

# ── Configuration ──────────────────────────────────────────────────────────────

SIGNIFICANCE_THRESHOLD = 1.0  # |log2FC| >= 1.0 = 2x fold-change

TIMEPOINT_WEIGHTS = {
    'R+1':   0.5,   # downweighted — reentry stress confound
    'R+45':  1.0,
    'R+82':  1.0,
    'R+194': 1.0,
}

# CBC metrics mapped to canonical biological domains.
# Only metrics with clear mechanistic relevance to each domain.
# Weights: 1.0 = primary indicator, 0.5 = secondary.
CBC_DOMAIN_MAP = {
    'immune_regulation': {
        'WHITE BLOOD CELL COUNT': 1.0,   # overall immune cell burden
        'LYMPHOCYTES':            1.0,   # adaptive immune arm
        'ABSOLUTE LYMPHOCYTES':   1.0,
        'EOSINOPHILS':            0.5,
        'ABSOLUTE EOSINOPHILS':   0.5,
        'BASOPHILS':              0.5,
        'ABSOLUTE BASOPHILS':     0.5,
    },
    'inflammation': {
        'NEUTROPHILS':            1.0,   # primary innate inflammatory cell
        'ABSOLUTE NEUTROPHILS':   1.0,
        'MONOCYTES':              0.5,   # secondary inflammatory mediator
        'ABSOLUTE MONOCYTES':     0.5,
        'PLATELET COUNT':         0.5,   # acute phase response marker
    },
    'oxidative_stress': {
        'RED BLOOD CELL COUNT':   1.0,   # space anemia — oxidative RBC destruction
        'HEMOGLOBIN':             1.0,
        'HEMATOCRIT':             1.0,
        'MCV':                    0.5,   # erythrocyte volume — erythropoietic stress proxy
        'RDW':                    0.5,   # red cell distribution width — anisocytosis
    },
}

# ── Helpers ────────────────────────────────────────────────────────────────────

def _fraction_to_risk(fraction):
    """Bin weighted fraction of significant markers to 1-5 risk score."""
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                  return 5


# ── Main ───────────────────────────────────────────────────────────────────────

def score_cbc(log2fc_df, subject_col='SUBJECT_ID', timepoint_col='TEST_DATE'):
    """
    Score CBC data per astronaut per domain.

    Parameters
    ----------
    log2fc_df    : DataFrame from CBC.py log2FC pipeline.
                   Columns: SUBJECT_ID, TEST_DATE, [CBC metric columns...].
                   Rows: one per astronaut x postflight timepoint.
    subject_col  : column name for subject ID
    timepoint_col: column name for timepoint

    Returns
    -------
    DataFrame with columns: astronaut, domain, overall_fraction, risk_score, source
    """
    records = []

    for domain, marker_weights in CBC_DOMAIN_MAP.items():
        # Only use markers that are actually in the dataframe
        available = {m: w for m, w in marker_weights.items() if m in log2fc_df.columns}
        if not available:
            print(f"  Warning: no markers found for domain '{domain}'")
            continue
        total_weight = sum(available.values())

        for astronaut in sorted(log2fc_df[subject_col].unique()):
            ast_data = log2fc_df[log2fc_df[subject_col] == astronaut]
            weighted_sum = 0.0
            weight_sum   = 0.0

            for _, row in ast_data.iterrows():
                tp   = row[timepoint_col]
                tp_w = TIMEPOINT_WEIGHTS.get(tp, 1.0)

                sig_w = sum(
                    w for m, w in available.items()
                    if pd.notna(row.get(m)) and abs(row.get(m, 0)) >= SIGNIFICANCE_THRESHOLD
                )
                fraction     = sig_w / total_weight if total_weight > 0 else 0.0
                weighted_sum += fraction * tp_w
                weight_sum   += tp_w

            overall_fraction = weighted_sum / weight_sum if weight_sum > 0 else 0.0

            records.append({
                'astronaut':        astronaut,
                'domain':           domain,
                'overall_fraction': round(overall_fraction, 4),
                'risk_score':       _fraction_to_risk(overall_fraction),
                'source':           'CBC',
            })

    df = pd.DataFrame(records)
    print(f"[CBC] Scored {len(df)} astronaut-domain pairs across {len(CBC_DOMAIN_MAP)} domains.")
    return df
