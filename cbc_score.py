"""
cbc_score.py — CBC domain scoring for AstronautMD
Scores using clinical reference ranges on raw values, not log2FC.
A marker is flagged if its raw value falls outside RANGE_MIN/RANGE_MAX
at a given postflight timepoint. Fraction of flagged markers = domain score.
"""

import numpy as np
import pandas as pd

TIMEPOINT_WEIGHTS = {
    'R+1':   0.5,
    'R+45':  1.0,
    'R+82':  1.0,
    'R+194': 1.0,
}

CBC_DOMAIN_MAP = {
    'immune_regulation': {
        'WHITE BLOOD CELL COUNT': 1.0,
        'LYMPHOCYTES':            1.0,
        'ABSOLUTE LYMPHOCYTES':   1.0,
        'EOSINOPHILS':            0.5,
        'ABSOLUTE EOSINOPHILS':   0.5,
        'BASOPHILS':              0.5,
        'ABSOLUTE BASOPHILS':     0.5,
    },
    'inflammation': {
        'NEUTROPHILS':            1.0,
        'ABSOLUTE NEUTROPHILS':   1.0,
        'MONOCYTES':              0.5,
        'ABSOLUTE MONOCYTES':     0.5,
        'PLATELET COUNT':         0.5,
    },
    'oxidative_stress': {
        'RED BLOOD CELL COUNT':   1.0,
        'HEMOGLOBIN':             1.0,
        'HEMATOCRIT':             1.0,
        'MCV':                    0.5,
        'RDW':                    0.5,
    },
}

def _fraction_to_risk(fraction):
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                  return 5


def score_cbc(cbc_raw, subject_col='SUBJECT_ID', timepoint_col='TEST_DATE'):
    """
    Score CBC data per astronaut per domain using clinical reference ranges.

    Parameters
    ----------
    cbc_raw      : raw CBC dataframe (long format) with columns:
                   ANALYTE, VALUE, RANGE_MIN, RANGE_MAX, SUBJECT_ID, TEST_DATE
    subject_col  : column name for subject ID
    timepoint_col: column name for timepoint

    Returns
    -------
    DataFrame with columns: astronaut, domain, overall_fraction, risk_score, source
    """
    postflight = ['R+1', 'R+45', 'R+82', 'R+194']

    # Build a lookup: (subject, timepoint, analyte) -> (value, range_min, range_max)
    post = cbc_raw[cbc_raw[timepoint_col].isin(postflight)].copy()
    post = post.reset_index(drop=True)

    # Flag out-of-range values
    post['out_of_range'] = (
        (post['VALUE'] < post['RANGE_MIN']) |
        (post['VALUE'] > post['RANGE_MAX'])
    )

    records = []

    for domain, marker_weights in CBC_DOMAIN_MAP.items():
        for astronaut in sorted(post[subject_col].unique()):
            ast = post[post[subject_col] == astronaut]
            weighted_sum = 0.0
            weight_sum   = 0.0

            for tp in postflight:
                tp_data = ast[ast[timepoint_col] == tp]
                if tp_data.empty:
                    continue

                tp_w         = TIMEPOINT_WEIGHTS.get(tp, 1.0)
                total_weight = 0.0
                flagged_weight = 0.0

                for marker, w in marker_weights.items():
                    row = tp_data[tp_data['ANALYTE'] == marker]
                    if row.empty:
                        continue
                    total_weight += w
                    if row['out_of_range'].values[0]:
                        flagged_weight += w

                if total_weight > 0:
                    fraction      = flagged_weight / total_weight
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