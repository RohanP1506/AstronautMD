"""
m6a_score.py — Per-astronaut m6A modification scoring for AstronautMD

Uses per-astronaut modification probability columns in m6a_filtered to
compute how many significant positions each crew member shows activity at
R+1 vs their own pre-flight baseline.

Biological rationale: m6A RNA methylation is a stress-responsive epigenetic
mark. Positions with increased modification at R+1 suggest post-transcriptional
regulation of stress response genes, mapped here to DNA damage response domain
(closest canonical domain given YTHDF reader roles in DNA damage signaling).
"""

import numpy as np
import pandas as pd

# ── Configuration ──────────────────────────────────────────────────────────────

PREFLIGHT_TIMEPOINTS  = ['L-92', 'L-44', 'L-3']
MOD_DIFF_THRESHOLD    = 0.10   # ≥10 percentage point shift in mod probability
CANONICAL_DOMAIN      = 'dna_damage'  # m6A contributes to DNA damage domain

# ── Helpers ────────────────────────────────────────────────────────────────────

def _fraction_to_risk(fraction):
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                  return 5


# ── Main ───────────────────────────────────────────────────────────────────────

def score_m6a(m6a_filtered):
    """
    Per-astronaut m6A activity score at R+1 across significant positions.

    Parameters
    ----------
    m6a_filtered : DataFrame of 94 significant m6A positions (from m6A.ipynb).
                   Must contain per-astronaut columns: C001_L-92, C001_R+1, etc.

    Returns
    -------
    DataFrame with columns: astronaut, domain, overall_fraction, risk_score,
                            source, n_sig_positions, n_total_valid_positions
    """
    subjects = ['C001', 'C002', 'C003', 'C004']
    records  = []

    for subject in subjects:
        pre_cols  = [f'{subject}_{tp}' for tp in PREFLIGHT_TIMEPOINTS
                     if f'{subject}_{tp}' in m6a_filtered.columns]
        post_col  = f'{subject}_R+1'

        if not pre_cols or post_col not in m6a_filtered.columns:
            print(f"  Warning: missing m6A columns for {subject}, skipping.")
            continue

        pre_mean = m6a_filtered[pre_cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)
        post_val = pd.to_numeric(m6a_filtered[post_col], errors='coerce')

        # Only positions where both pre and post are available
        valid_mask = pre_mean.notna() & post_val.notna()
        n_valid    = valid_mask.sum()

        if n_valid == 0:
            print(f"  Warning: no valid positions for {subject}, skipping.")
            continue

        diff  = (post_val[valid_mask] - pre_mean[valid_mask]).abs()
        n_sig = int((diff >= MOD_DIFF_THRESHOLD).sum())

        fraction = n_sig / n_valid if n_valid > 0 else 0.0

        records.append({
            'astronaut':               subject,
            'domain':                  CANONICAL_DOMAIN,
            'overall_fraction':        round(fraction, 4),
            'risk_score':              _fraction_to_risk(fraction),
            'source':                  'm6A',
            'n_sig_positions':         n_sig,
            'n_total_valid_positions': int(n_valid),
        })

    df = pd.DataFrame(records)
    print(f"[m6A] Scored {len(df)} astronauts on domain '{CANONICAL_DOMAIN}'.")
    return df
