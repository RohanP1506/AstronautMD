"""
rnaseq_score.py — Per-astronaut RNA-seq scoring for AstronautMD

DEGs in degs_FP1.csv are cohort-level (DESeq2 pooled across all 4 astronauts),
but the file also contains raw featureCounts per astronaut per timepoint.
This module uses those raw counts to compute a per-astronaut log2FC for each
DEG at R+1 vs that astronaut's own pre-flight mean, then scores by domain.

Important caveat: single-sample per timepoint — no statistical testing per
astronaut. These are magnitude estimates only, not significance calls.
"""

import numpy as np
import pandas as pd

# ── Configuration ──────────────────────────────────────────────────────────────

SIGNIFICANCE_THRESHOLD = 1.0  # |log2FC| >= 1.0 for per-astronaut magnitude
PREFLIGHT_TIMEPOINTS   = ['L-92', 'L-44', 'L-3']
TIMEPOINT_WEIGHT_R1    = 0.5   # R+1 is only available timepoint, downweighted

# Map RNA-seq domain labels (from degs_FP1.csv) to canonical 5-domain system.
# 'Other' is excluded — too heterogeneous to assign meaningfully.
RNASEQ_TO_CANONICAL = {
    'Erythroid_RBC':              'oxidative_stress',    # space anemia = oxidative RBC destruction
    'Adaptive_Immunity':          'immune_regulation',
    'Mitochondrial':              'mitochondrial',
    'Protein_Degradation':        'dna_damage',          # ubiquitin/proteasome = damage response
    'Transcriptional_Regulation': 'immune_regulation',   # transcription factors active in immune cells
    'Other':                      None,                  # skip
}

# ── Helpers ────────────────────────────────────────────────────────────────────

def _fraction_to_risk(fraction):
    if fraction < 0.10:   return 1
    elif fraction < 0.25: return 2
    elif fraction < 0.50: return 3
    elif fraction < 0.75: return 4
    else:                  return 5


# ── Main ───────────────────────────────────────────────────────────────────────

def score_rnaseq(degs_path='observations/whole_blood/degs_FP1.csv'):
    """
    Compute per-astronaut RNA-seq domain scores from featureCounts.

    Parameters
    ----------
    degs_path : path to degs_FP1.csv (exported from RNAseq.ipynb)

    Returns
    -------
    DataFrame with columns: astronaut, domain, overall_fraction, risk_score,
                            source, rna_domain, n_genes
    """
    degs = pd.read_csv(degs_path, index_col=0)

    if 'domain' not in degs.columns:
        raise ValueError("degs_FP1.csv must have a 'domain' column. "
                         "Run domain mapping step in RNAseq.ipynb first.")

    subjects = ['C001', 'C002', 'C003', 'C004']
    records  = []

    for subject in subjects:
        # Build pre-flight and post-flight column names
        pre_cols  = [f'{subject}_{tp}' for tp in PREFLIGHT_TIMEPOINTS
                     if f'{subject}_{tp}' in degs.columns]
        post_col  = f'{subject}_R+1'

        if not pre_cols or post_col not in degs.columns:
            print(f"  Warning: missing featureCounts columns for {subject}, skipping.")
            continue

        # Per-gene log2FC for this astronaut
        pre_mean = degs[pre_cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)
        post_val = pd.to_numeric(degs[post_col], errors='coerce')

        valid_mask = (pre_mean > 0) & (post_val > 0)
        gene_log2fc = pd.Series(np.nan, index=degs.index)
        gene_log2fc[valid_mask] = np.log2(post_val[valid_mask] / pre_mean[valid_mask])

        # Score per RNA-seq domain
        for rna_domain in degs['domain'].unique():
            canonical = RNASEQ_TO_CANONICAL.get(rna_domain)
            if canonical is None:
                continue

            domain_genes = degs[degs['domain'] == rna_domain]
            domain_fc    = gene_log2fc[domain_genes.index].dropna()

            if len(domain_fc) == 0:
                continue

            sig_frac         = (domain_fc.abs() >= SIGNIFICANCE_THRESHOLD).sum() / len(domain_fc)
            overall_fraction = float(sig_frac) * TIMEPOINT_WEIGHT_R1

            records.append({
                'astronaut':        subject,
                'domain':           canonical,
                'overall_fraction': round(overall_fraction, 4),
                'risk_score':       _fraction_to_risk(overall_fraction),
                'source':           'RNAseq',
                'rna_domain':       rna_domain,
                'n_genes':          len(domain_fc),
            })

    df = pd.DataFrame(records)
    print(f"[RNA-seq] Scored {len(df)} astronaut-domain pairs.")
    return df
