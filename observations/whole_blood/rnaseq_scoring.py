import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# WHOLE BLOOD RNA-seq DOMAIN SCORING
# OSD-569 / GLDS-561  —  only R+1 timepoint available
#
# Methodology (different from cytokine scoring — no trajectory possible):
#   - Filter DEGs: padj < 0.05 AND |log2FC| >= 1.0 using DESeq2 results
#   - Assign each DEG to a biological domain via keyword matching on gene sets
#   - Domain score = weighted sum of |log2FC| for significant genes in domain,
#     normalized by total possible weight in that domain
#   - Risk tier: normalized score binned to 1-5 scale
#   - Note: R+1 only — this is a snapshot, not a trajectory
#
# Why not the same scoring as cytokines?
#   Cytokine scoring uses fraction of markers significant (bounded 0-1).
#   RNAseq has hundreds of genes per domain — fraction would be near-zero.
#   Instead we use weighted mean |log2FC| of significant DEGs, which captures
#   both the number of genes affected and the magnitude of their dysregulation.
# =============================================================================

# Significance thresholds
PADJ_THRESHOLD  = 0.05
LOG2FC_THRESHOLD = 1.0   # 2x fold-change

# =============================================================================
# DOMAIN GENE SETS
# Keywords matched against Ensembl gene IDs or gene symbols.
# Based on g:Profiler enrichment results from RNAseq notebook.
# Domains mirror the ones used for the overall astronaut risk profile.
# =============================================================================

RNASEQ_DOMAINS = {
    'Anemia_Risk': {
        'genes': [
            # Erythroid structural and regulatory genes from DEG analysis
            'HBB', 'HBA1', 'HBA2', 'GYPA', 'GYPB', 'GYPE',   # hemoglobin + glycophorins
            'ALAS2', 'SLC4A1', 'EPB42', 'EKLF', 'KLF1',       # erythropoiesis regulators
            'TFRC', 'SLC25A37', 'FECH', 'HMBS',                # iron/heme metabolism
            'RHAG', 'RHCE', 'RHD', 'ANK1', 'SPTA1', 'SPTB',   # RBC membrane
        ],
        'weight': 1.0,
        'description': 'Red blood cell biology and space anemia risk'
    },
    'Immune_Dysregulation': {
        'genes': [
            # Adaptive immunity genes
            'CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B',
            'GZMB', 'GZMA', 'PRF1',                            # cytotoxic T cell
            'CD19', 'CD79A', 'CD79B',                           # B cell
            'FOXP3', 'IL2RA',                                   # Treg
            'IFNG', 'TNF', 'IL2',                               # effector cytokines
            'KLRD1', 'KLRK1', 'NKG7',                          # NK cell
        ],
        'weight': 1.0,
        'description': 'Adaptive immune system dysregulation'
    },
    'Mitochondrial_Stress': {
        'genes': [
            # Mitochondrial function genes
            'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ND1', 'MT-ND2',
            'MT-ATP6', 'MT-CYB',
            'UQCRQ', 'UQCRC1', 'UQCRC2',
            'COX4I1', 'COX5A', 'COX6A1',
            'NDUFA1', 'NDUFB3', 'NDUFC1',
            'VDAC1', 'VDAC2', 'TOMM20',
        ],
        'weight': 0.8,
        'description': 'Mitochondrial function and cellular energy stress'
    },
    'Proteostasis': {
        'genes': [
            # Ubiquitin-proteasome system
            'UBB', 'UBC', 'UBA52',
            'PSMA1', 'PSMA2', 'PSMB1', 'PSMB2', 'PSMB3',
            'PSMC1', 'PSMC2', 'PSMD1', 'PSMD2',
            'UCHL1', 'UCHL3',
            'USP14', 'USP7',
            'HSPA1A', 'HSPA1B', 'HSP90AA1', 'DNAJB1',         # chaperones
        ],
        'weight': 0.8,
        'description': 'Protein quality control and degradation stress'
    },
    'Transcriptional_Stress': {
        'genes': [
            # Transcription factors and chromatin regulators
            'EGR1', 'EGR2', 'EGR3',
            'FOS', 'FOSB', 'JUN', 'JUNB',
            'ATF3', 'ATF4',
            'KLF4', 'KLF6',
            'NR4A1', 'NR4A2', 'NR4A3',
            'STAT1', 'STAT3',
            'MYC', 'MYCN',
        ],
        'weight': 0.6,
        'description': 'Stress-responsive transcriptional reprogramming'
    },
}


# =============================================================================
# DATA LOADING
# =============================================================================

def load_rnaseq(url=None, filepath=None):
    """
    Load RNAseq data from URL or local file.
    Expects the GLDS-561 format with DESeq2 columns.
    """
    if url:
        df = pd.read_excel(url, index_col=0)
    elif filepath:
        if filepath.endswith('.csv'):
            df = pd.read_csv(filepath, index_col=0)
        else:
            df = pd.read_excel(filepath, index_col=0)
    else:
        raise ValueError("Provide either url or filepath")
    return df


def extract_gene_symbol(ensembl_id):
    """
    Placeholder — in practice you'd use a gene annotation table.
    Returns the Ensembl ID stripped of version number.
    """
    return ensembl_id.split('.')[0]


# =============================================================================
# SCORING
# =============================================================================

def score_rnaseq_domains(rnaseq_df, gene_symbol_map=None):
    """
    Score each RNAseq domain per astronaut at R+1.

    Parameters
    ----------
    rnaseq_df       : raw RNAseq dataframe (GLDS-561 format)
    gene_symbol_map : dict mapping Ensembl ID -> gene symbol (optional).
                      If None, uses Ensembl IDs for matching.

    Returns
    -------
    domain_scores   : DataFrame with columns [astronaut, domain, risk_score,
                      n_sig_genes, mean_abs_log2fc, description]
    sig_genes_df    : DataFrame of all significant DEGs with domain assignment
    """
    # Extract per-astronaut R+1 log2FC
    # Columns are like C001_R+1, C001_L-92 etc.
    astronauts = ['C001', 'C002', 'C003', 'C004']
    records = []

    for astro in astronauts:
        # Build personal preflight baseline (mean across L- timepoints)
        preflight_cols = [c for c in rnaseq_df.columns if f'{astro}_L-' in c and '.1' not in c]
        postflight_col = f'{astro}_R+1'

        if postflight_col not in rnaseq_df.columns or not preflight_cols:
            continue

        pre_mean  = rnaseq_df[preflight_cols].mean(axis=1)
        post_vals = rnaseq_df[postflight_col]

        # Compute personal log2FC (guard against zeros)
        valid = (pre_mean > 0) & (post_vals > 0)
        log2fc = pd.Series(np.nan, index=rnaseq_df.index)
        log2fc[valid] = np.log2(post_vals[valid] / pre_mean[valid])

        # Also use DESeq2 padj for population-level significance filter
        padj_col = 'DESeq2_adjusted p-value'
        padj = rnaseq_df[padj_col] if padj_col in rnaseq_df.columns else pd.Series(1.0, index=rnaseq_df.index)

        for gene_id in rnaseq_df.index:
            symbol = gene_symbol_map.get(gene_id, gene_id) if gene_symbol_map else gene_id
            gene_log2fc = log2fc.get(gene_id, np.nan)
            gene_padj   = padj.get(gene_id, 1.0)

            records.append({
                'astronaut':   astro,
                'gene_id':     gene_id,
                'gene_symbol': symbol,
                'log2fc':      gene_log2fc,
                'padj':        gene_padj,
            })

    gene_df = pd.DataFrame(records)

    # Flag significance: use DESeq2 padj for reliability
    gene_df['significant'] = (
        (gene_df['padj'] < PADJ_THRESHOLD) &
        (gene_df['log2fc'].abs() >= LOG2FC_THRESHOLD)
    )

    # Assign domains by gene symbol matching
    def assign_domain(symbol):
        symbol_upper = str(symbol).upper()
        for domain, info in RNASEQ_DOMAINS.items():
            if any(g.upper() == symbol_upper or g.upper() in symbol_upper
                   for g in info['genes']):
                return domain
        return None

    gene_df['domain'] = gene_df['gene_symbol'].apply(assign_domain)
    sig_genes_df = gene_df[gene_df['significant'] & gene_df['domain'].notna()].copy()

    # Score per astronaut per domain
    domain_records = []
    for astro in astronauts:
        for domain, info in RNASEQ_DOMAINS.items():
            subset = sig_genes_df[
                (sig_genes_df['astronaut'] == astro) &
                (sig_genes_df['domain'] == domain)
            ]

            n_sig = len(subset)
            mean_abs_log2fc = subset['log2fc'].abs().mean() if n_sig > 0 else 0.0

            # Risk score: based on mean |log2FC| of significant domain genes
            # 0 genes = 1, small effect = 2, moderate = 3, large = 4, severe = 5
            if n_sig == 0:
                risk = 1
            elif mean_abs_log2fc < 1.5:
                risk = 2
            elif mean_abs_log2fc < 2.0:
                risk = 3
            elif mean_abs_log2fc < 3.0:
                risk = 4
            else:
                risk = 5

            domain_records.append({
                'astronaut':       astro,
                'domain':          domain,
                'risk_score':      risk,
                'n_sig_genes':     n_sig,
                'mean_abs_log2fc': round(mean_abs_log2fc, 3),
                'description':     info['description'],
            })

    domain_scores = pd.DataFrame(domain_records)

    # Printout
    print(f"\n{'='*65}")
    print(f"  Whole Blood RNA-seq Domain Scores (R+1, OSD-569)")
    print(f"  Threshold: padj < {PADJ_THRESHOLD}, |log2FC| >= {LOG2FC_THRESHOLD}")
    print(f"  Note: R+1 only — snapshot score, no trajectory available")
    print(f"{'='*65}")

    for domain in RNASEQ_DOMAINS:
        print(f"\n  {domain.replace('_', ' ')}")
        subset = domain_scores[domain_scores['domain'] == domain]
        print(subset[['astronaut', 'n_sig_genes', 'mean_abs_log2fc', 'risk_score']].to_string(index=False))

    return domain_scores, sig_genes_df


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_rnaseq_domain_heatmap(domain_scores, output_path='fig_rnaseq_domain_scores.png'):
    """
    Heatmap of RNAseq risk scores: rows = domains, columns = astronauts.
    """
    pivot = domain_scores.pivot_table(
        index='domain', columns='astronaut', values='risk_score'
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.heatmap(
        pivot, ax=ax,
        cmap='YlOrRd', vmin=1, vmax=5,
        annot=True, fmt='.0f', annot_kws={'size': 11, 'weight': 'bold'},
        linewidths=0.5, linecolor='#dddddd',
        cbar_kws={'label': 'Risk Score (1–5)', 'shrink': 0.7}
    )
    ax.set_title(
        'Whole Blood RNA-seq Domain Risk Scores\n'
        'OSD-569 — Day R+1 vs Preflight Baseline',
        fontsize=12, fontweight='bold'
    )
    ax.set_xlabel('Astronaut')
    ax.set_ylabel('')
    ax.set_yticklabels(
        [d.replace('_', ' ') for d in pivot.index],
        rotation=0, fontsize=9
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f'[✓] Saved: {output_path}')
    plt.close(fig)


def plot_rnaseq_deg_bar(sig_genes_df, output_path='fig_rnaseq_deg_counts.png'):
    """
    Horizontal bar chart: number of significant DEGs per domain per astronaut.
    """
    counts = (
        sig_genes_df
        .groupby(['domain', 'astronaut'])
        .size()
        .reset_index(name='n_degs')
    )

    domains = list(RNASEQ_DOMAINS.keys())
    astronauts = sorted(sig_genes_df['astronaut'].unique())
    colors = plt.cm.tab10.colors

    fig, axes = plt.subplots(1, len(domains), figsize=(4 * len(domains), 4), sharey=False)
    if len(domains) == 1:
        axes = [axes]

    for ax, domain in zip(axes, domains):
        sub = counts[counts['domain'] == domain]
        vals = [sub.loc[sub['astronaut'] == a, 'n_degs'].values[0]
                if a in sub['astronaut'].values else 0
                for a in astronauts]
        bars = ax.bar(astronauts, vals,
                      color=[colors[i % len(colors)] for i in range(len(astronauts))],
                      edgecolor='white', linewidth=0.5)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                    str(val), ha='center', fontsize=9, fontweight='bold')
        ax.set_title(domain.replace('_', ' '), fontsize=9, fontweight='bold')
        ax.set_xlabel('Astronaut')
        ax.set_ylabel('Significant DEGs' if ax == axes[0] else '')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(
        'Significant DEGs per Domain at R+1\nWhole Blood RNA-seq (OSD-569)',
        fontsize=12, fontweight='bold'
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f'[✓] Saved: {output_path}')
    plt.close(fig)