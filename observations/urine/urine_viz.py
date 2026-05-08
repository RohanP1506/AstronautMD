import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import re

# ─────────────────────────────────────────────
# TIMEPOINT SORTING UTILITY
# Converts L-180 → -180, L-90 → -90, R+0 → 0, R+30 → 30
# so timepoints sort chronologically on any axis.
# ─────────────────────────────────────────────
def timepoint_to_int(tp):
    tp = str(tp).strip()
    m = re.match(r'L-(\d+)', tp)
    if m:
        return -int(m.group(1))
    m = re.match(r'R\+(\d+)', tp)
    if m:
        return int(m.group(1))
    return 0  # fallback (e.g. 'FD' or unknown)

def sort_timepoints(timepoints):
    return sorted(timepoints, key=timepoint_to_int)

# ─────────────────────────────────────────────
# STYLE SETUP
# ─────────────────────────────────────────────
PALETTE = {
    'inflammation':    '#E05C5C',
    'oxidative_stress':'#5B8DD9',
}

ASTRONAUT_COLORS = sns.color_palette('tab10')

plt.rcParams.update({
    'font.family':      'DejaVu Sans',
    'axes.spines.top':  False,
    'axes.spines.right':False,
    'axes.titlesize':   11,
    'axes.labelsize':   10,
    'xtick.labelsize':  8,
    'ytick.labelsize':  8,
    'figure.dpi':       150,
})

# ─────────────────────────────────────────────
# PLOT 1: COMPOSITE SCORE TIMELINE
# One panel per domain; lines = astronauts
# ─────────────────────────────────────────────
def plot_composite_timeline(composite_scores, output_path='fig1_composite_timeline.png'):
    """
    Line plot: composite_score over chronological timepoints.
    Faceted by domain. Each line = one astronaut.
    A vertical dashed line marks launch (between last L- and R+0/R+1).
    """
    df = composite_scores.copy()
    df['tp_int'] = df['timepoint'].apply(timepoint_to_int)

    all_timepoints = sort_timepoints(df['timepoint'].unique())
    tp_order = {tp: i for i, tp in enumerate(all_timepoints)}
    df['tp_order'] = df['timepoint'].map(tp_order)

    domains = df['domain'].unique()
    astronauts = sorted(df['astronaut'].unique())
    n_domains = len(domains)

    fig, axes = plt.subplots(1, n_domains, figsize=(6 * n_domains, 4.5), sharey=True)
    if n_domains == 1:
        axes = [axes]

    # Find launch line position (between last negative and first non-negative)
    tp_ints = sorted(df['tp_int'].unique())
    launch_x = None
    for i in range(len(all_timepoints) - 1):
        a = timepoint_to_int(all_timepoints[i])
        b = timepoint_to_int(all_timepoints[i + 1])
        if a < 0 and b >= 0:
            launch_x = i + 0.5
            break

    for ax, domain in zip(axes, sorted(domains)):
        sub = df[df['domain'] == domain]
        color = PALETTE.get(domain, '#888888')

        for j, astro in enumerate(astronauts):
            asub = sub[sub['astronaut'] == astro].sort_values('tp_order')
            if asub.empty:
                continue
            ax.plot(
                asub['tp_order'], asub['composite_score'],
                marker='o', markersize=5, linewidth=1.8,
                color=ASTRONAUT_COLORS[j % len(ASTRONAUT_COLORS)],
                label=astro, alpha=0.85
            )

        if launch_x is not None:
            ax.axvline(launch_x, color='black', linestyle='--', linewidth=1.2, alpha=0.5, label='Launch')

        ax.set_title(domain.replace('_', ' ').title(), color=color, fontweight='bold')
        ax.set_xticks(range(len(all_timepoints)))
        ax.set_xticklabels(all_timepoints, rotation=45, ha='right')
        ax.set_ylim(0, 5.2)
        ax.set_xlabel('Timepoint')
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.axhspan(2.0, 5.2, color=color, alpha=0.05)  # highlight elevated zone

    axes[0].set_ylabel('Composite Score (0–5)')
    axes[-1].legend(title='Astronaut', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)

    fig.suptitle('Urine Immune Domain Composite Scores\nOSD-656 Alamar Panel', 
                 fontsize=13, fontweight='bold', y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    print(f"[✓] Saved: {output_path}")
    plt.close(fig)


# ─────────────────────────────────────────────
# PLOT 2: LOG2FC HEATMAP PER ASTRONAUT
# Rows = markers, Columns = timepoints
# Diverging color: blue=down, red=up
# ─────────────────────────────────────────────
def plot_log2fc_heatmap(urine_processed, output_path='fig2_log2fc_heatmap.png',
                        min_abs_fc=0.5):
    """
    Heatmap of log2FC. One subplot per astronaut.
    Only shows markers that hit |log2FC| >= min_abs_fc in at least one timepoint.
    """
    df = urine_processed.copy()

    # Clean marker names for display (strip long suffixes like '_concentration_NPQ')
    df['marker_short'] = df['marker'].str.extract(r'^([^_]+(?:_[^_]+)?)')[0].str.upper()

    # Keep only markers with meaningful fold change somewhere
    max_abs = df.groupby('marker_short')['log2fc'].apply(lambda x: x.abs().max())
    keep = max_abs[max_abs >= min_abs_fc].index
    df = df[df['marker_short'].isin(keep)]

    astronauts = sorted(df['astronaut'].unique())
    all_timepoints = sort_timepoints(df['timepoint'].unique())
    all_markers = sorted(df['marker_short'].unique())

    n = len(astronauts)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, max(4, 0.35 * len(all_markers) + 1.5)),
                             sharey=True)
    if n == 1:
        axes = [axes]

    vmax = min(df['log2fc'].abs().quantile(0.95), 3.5)  # cap color scale at 95th percentile

    for ax, astro in zip(axes, astronauts):
        sub = df[df['astronaut'] == astro]
        pivot = sub.pivot_table(index='marker_short', columns='timepoint',
                                values='log2fc', aggfunc='mean')
        # Reindex columns to chronological order (only columns present in data)
        ordered_cols = [tp for tp in all_timepoints if tp in pivot.columns]
        pivot = pivot.reindex(columns=ordered_cols, index=all_markers)

        sns.heatmap(
            pivot, ax=ax,
            cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.3, linecolor='#dddddd',
            cbar=(ax == axes[-1]),
            cbar_kws={'label': 'log₂FC', 'shrink': 0.6} if ax == axes[-1] else {}
        )
        ax.set_title(astro, fontweight='bold')
        ax.set_xlabel('Timepoint')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=7)
        if ax == axes[0]:
            ax.set_ylabel('Marker')
        else:
            ax.set_ylabel('')

    fig.suptitle('Urine Cytokine Log₂FC Relative to Preflight Baseline\nOSD-656 Alamar Panel',
                 fontsize=13, fontweight='bold', y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    print(f"[✓] Saved: {output_path}")
    plt.close(fig)


# ─────────────────────────────────────────────
# PLOT 3: POST-FLIGHT PERSISTENCE BAR CHART
# How many R+ timepoints each marker stayed elevated
# Faceted per astronaut
# ─────────────────────────────────────────────
def plot_persistence_bars(urine_processed, output_path='fig3_persistence_bars.png',
                           top_n=10):
    """
    Horizontal bar chart: markers ranked by how many post-flight
    timepoints they had marker_score >= 1.
    Color = mean log2FC direction (red=up, blue=down).
    One subplot per astronaut.
    """
    df = urine_processed.copy()
    df['marker_short'] = df['marker'].str.extract(r'^([^_]+(?:_[^_]+)?)')[0].str.upper()

    postflight = df[df['timepoint'].str.startswith('R+')]
    sig = postflight[postflight['marker_score'] >= 1]

    persistence = sig.groupby(['astronaut', 'marker_short']).agg(
        times_significant=('marker_score', 'count'),
        mean_log2fc=('log2fc', 'mean'),
        max_score=('marker_score', 'max')
    ).reset_index()

    astronauts = sorted(persistence['astronaut'].unique())
    n = len(astronauts)

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, astro in zip(axes, astronauts):
        sub = persistence[persistence['astronaut'] == astro]
        sub = sub.sort_values(['times_significant', 'max_score'], ascending=True).tail(top_n)

        if sub.empty:
            ax.text(0.5, 0.5, 'No persistent markers', ha='center', va='center', 
                    transform=ax.transAxes, fontsize=10, color='gray')
            ax.set_title(astro, fontweight='bold')
            continue

        colors = ['#E05C5C' if v > 0 else '#5B8DD9' for v in sub['mean_log2fc']]
        bars = ax.barh(sub['marker_short'], sub['times_significant'],
                       color=colors, edgecolor='white', linewidth=0.5, height=0.65)

        # Annotate with mean log2FC
        for bar, (_, row) in zip(bars, sub.iterrows()):
            ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                    f"{row['mean_log2fc']:+.2f}", va='center', fontsize=7.5, color='#333333')

        ax.set_xlabel('# R+ Timepoints Significant')
        ax.set_title(astro, fontweight='bold')
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.set_xlim(0, sub['times_significant'].max() + 1.5)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E05C5C', label='Upregulated (mean log₂FC > 0)'),
        Patch(facecolor='#5B8DD9', label='Downregulated (mean log₂FC < 0)')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2,
               bbox_to_anchor=(0.5, -0.05), fontsize=9, frameon=False)

    fig.suptitle('Post-Flight Marker Persistence\nOSD-656 Urine Cytokines (R+ Timepoints, Score ≥ 1)',
                 fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    print(f"[✓] Saved: {output_path}")
    plt.close(fig)


# ─────────────────────────────────────────────
# MAIN — call this from urine.py
# ─────────────────────────────────────────────
def generate_all_figures(urine_processed, composite_scores, output_dir='.'):
    import os
    os.makedirs(output_dir, exist_ok=True)

    plot_composite_timeline(
        composite_scores,
        output_path=f'{output_dir}/fig1_composite_timeline.png'
    )
    plot_log2fc_heatmap(
        urine_processed,
        output_path=f'{output_dir}/fig2_log2fc_heatmap.png'
    )
    plot_persistence_bars(
        urine_processed,
        output_path=f'{output_dir}/fig3_persistence_bars.png'
    )
    print(f"\n[✓] All figures saved to: {output_dir}/")