"""
cytokine_score.py — Serum cytokine scoring for AstronautMD
Wraps the existing immune_eve pipeline and returns standardized score dataframes.
"""
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'observations', 'cytokine'))

from observations.cytokine.cytokine_dataframe import (
    organize_df, process_cytokine_df, score, set_marker_weights
)

CYTOKINE_URL = (
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download'
    '?source=datamanager&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv'
)

MARKERS = {
    'immune_regulation': [
        'il_2_concentration_picogram_per_milliliter',
        'il_4_concentration_picogram_per_milliliter',
        'il_10_concentration_picogram_per_milliliter',
        'il_13_concentration_picogram_per_milliliter',
        'il_15_concentration_picogram_per_milliliter',
        'il_17a_concentration_picogram_per_milliliter',
        'il_1ra_concentration_picogram_per_milliliter',
        'il_21_concentration_picogram_per_milliliter',
        'il_27_concentration_picogram_per_milliliter',
        'ifny_concentration_picogram_per_milliliter',
    ],
    'inflammation': [
        'il_1a_concentration_picogram_per_milliliter',
        'il_1b_concentration_picogram_per_milliliter',
        'il_6_concentration_picogram_per_milliliter',
        'il_8_concentration_picogram_per_milliliter',
        'il_18_concentration_picogram_per_milliliter',
        'il_33_concentration_picogram_per_milliliter',
        'tnfa_concentration_picogram_per_milliliter',
        'tnfb_concentration_picogram_per_milliliter',
        'groa_concentration_picogram_per_milliliter',
        'mcp_1_concentration_picogram_per_milliliter',
        'mip_1a_concentration_picogram_per_milliliter',
        'mip_1b_concentration_picogram_per_milliliter',
        'ip_10_concentration_picogram_per_milliliter',
    ],
}

MARKER_WEIGHTS = {
    'il_2_concentration_picogram_per_milliliter':   1.0,
    'il_4_concentration_picogram_per_milliliter':   0.5,
    'il_10_concentration_picogram_per_milliliter':  1.0,
    'il_13_concentration_picogram_per_milliliter':  0.5,
    'il_15_concentration_picogram_per_milliliter':  1.0,
    'il_17a_concentration_picogram_per_milliliter': 0.5,
    'il_1ra_concentration_picogram_per_milliliter': 1.0,
    'il_21_concentration_picogram_per_milliliter':  0.5,
    'il_27_concentration_picogram_per_milliliter':  0.5,
    'ifny_concentration_picogram_per_milliliter':   1.0,
    'il_1a_concentration_picogram_per_milliliter':  0.5,
    'il_1b_concentration_picogram_per_milliliter':  1.0,
    'il_6_concentration_picogram_per_milliliter':   1.0,
    'il_8_concentration_picogram_per_milliliter':   1.0,
    'il_18_concentration_picogram_per_milliliter':  0.5,
    'il_33_concentration_picogram_per_milliliter':  0.5,
    'tnfa_concentration_picogram_per_milliliter':   1.0,
    'tnfb_concentration_picogram_per_milliliter':   0.5,
    'groa_concentration_picogram_per_milliliter':   0.5,
    'mcp_1_concentration_picogram_per_milliliter':  1.0,
    'mip_1a_concentration_picogram_per_milliliter': 0.5,
    'mip_1b_concentration_picogram_per_milliliter': 0.5,
    'ip_10_concentration_picogram_per_milliliter':  1.0,
}


def score_cytokines():
    """
    Load serum cytokine data, score immune_regulation and inflammation domains.
    Returns standardized DataFrame: astronaut, domain, overall_fraction,
                                    risk_score, source
    """
    raw = pd.read_csv(CYTOKINE_URL, index_col=0)
    long = organize_df(raw.transpose())
    processed = process_cytokine_df(long, 'immune_eve')
    set_marker_weights(MARKER_WEIGHTS)

    records = []
    for domain, domain_markers in MARKERS.items():
        _, overall = score(processed, domain_markers, domain_name=domain)
        if overall.empty:
            continue
        overall = overall.copy()
        overall['domain'] = domain
        overall['source'] = 'Cytokine'
        records.append(overall[['astronaut', 'domain', 'overall_fraction',
                                 'risk_score', 'source']])

    if not records:
        print("  WARNING: No cytokine scores produced.")
        return pd.DataFrame()

    df = pd.concat(records, ignore_index=True)
    print(f"[Cytokine] Scored {len(df)} astronaut-domain pairs.")
    return df