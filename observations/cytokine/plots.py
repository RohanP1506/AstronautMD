from immune_eve import *
import matplotlib.pyplot as plt
import seaborn as sns
import re


def sort_timepoints(timepoints):
    def key(tp):
        m = re.match(r'(L|R)[+-](\d+)', tp)
        if not m:
            return (1, 0)
        direction, num = m.group(1), int(m.group(2))
        return (0, -num) if direction == 'L' else (1, num)
    return sorted(timepoints, key=key)


def plot_log2fc_heatmap(processed_df, title_suffix='OSD-575 EVE Panel',
                        output_path='fig_log2fc_heatmap.png', min_abs_fc=0.5):
    df = processed_df.copy()

    # FIXED: full name stripping instead of token extraction — avoids collisions
    df['marker_short'] = (
        df['marker']
        .str.replace('_concentration_picogram_per_milliliter', '', regex=False)
        .str.upper()
    )

    max_abs = df.groupby('marker_short')['log2fc'].apply(lambda x: x.abs().max())
    keep = max_abs[max_abs >= min_abs_fc].index
    df = df[df['marker_short'].isin(keep)]

    astronauts = sorted(df['astronaut'].unique())
    all_timepoints = sort_timepoints(df['timepoint'].unique())
    all_markers = sorted(df['marker_short'].unique())

    n = len(astronauts)
    fig, axes = plt.subplots(
        1, n,
        figsize=(4.5 * n, max(4, 0.35 * len(all_markers) + 1.5)),
        sharey=True
    )
    if n == 1:
        axes = [axes]

    vmax = min(df['log2fc'].abs().quantile(0.95), 3.5)

    for ax, astro in zip(axes, astronauts):
        sub = df[df['astronaut'] == astro]
        pivot = sub.pivot_table(
            index='marker_short',
            columns='timepoint',
            values='log2fc',
            aggfunc='mean'
        )
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
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index, rotation=0, fontsize=7)
        ax.set_ylabel('Marker' if ax == axes[0] else '')

    fig.suptitle(
        f'Cytokine Log₂FC Relative to Preflight Baseline\n{title_suffix}',
        fontsize=13, fontweight='bold', y=1.05
    )
    fig.tight_layout()
    fig.subplots_adjust(top=0.82)
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f'[✓] Saved: {output_path}')
    plt.close(fig)


# ===================
# Log2FC heatmaps (one per domain)
# ===================

for domain_name, domain_markers in markers.items():
    domain_df = immune_eve_long[immune_eve_long['marker'].isin(domain_markers)]
    plot_log2fc_heatmap(
        domain_df,
        title_suffix=f'OSD-575 EVE Panel — {domain_name.replace("_", " ").title()}',
        output_path=f'fig_log2fc_heatmap_{domain_name}.png'
    )


# ===================
# Score vs. Time plot
# ===================

# for composite, long_df, label in [
#     (immune_reg_composite,    immune_eve_long, 'Immune Regulation'),
#     (inflammation_composite,  immune_eve_long, 'Inflammation')
# ]:
#     composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
#     composite = composite.sort_values(['astronaut', 'timepoint_num'])
#     timepoint_order = (
#         composite.drop_duplicates('timepoint')
#         .sort_values('timepoint_num')['timepoint']
#         .tolist()
#     )

#     # Build matching log2fc composite for this domain's markers
#     domain_markers = (
#         markers['immune_regulation'] if label == 'Immune Regulation'
#         else markers['inflammation']
#     )
#     log2fc_composite = (
#         long_df[
#             long_df['timepoint'].str.startswith('R+') &
#             long_df['marker'].isin(domain_markers)
#         ]
#         .groupby(['astronaut', 'timepoint'])['log2fc']
#         .mean()
#         .reset_index()
#     )
#     log2fc_composite['timepoint_num'] = (
#         log2fc_composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
#     )
#     log2fc_composite = log2fc_composite.sort_values(['astronaut', 'timepoint_num'])

#     fig, (ax_score, ax_log2fc) = plt.subplots(
#         2, 1, figsize=(11, 9), sharex=True
#     )

#     for astronaut, data in composite.groupby('astronaut'):
#         data = data.set_index('timepoint').reindex(timepoint_order)
#         ax_score.plot(timepoint_order, data['marker_score'],
#                       marker='o', linewidth=2, markersize=6, label=astronaut)

#     for astronaut, data in log2fc_composite.groupby('astronaut'):
#         data = data.set_index('timepoint').reindex(timepoint_order)
#         ax_log2fc.plot(timepoint_order, data['log2fc'],
#                        marker='o', linewidth=2, markersize=6, label=astronaut)

#     # --- marker score panel ---
#     ax_score.axhline(0, color='black', linestyle='--', linewidth=1.2,
#                      alpha=0.6, label='Baseline')
#     ax_score.axhspan(0, 5,  alpha=0.04, color='red')
#     ax_score.axhspan(-5, 0, alpha=0.04, color='blue')
#     ax_score.set_ylabel('Composite marker score', fontsize=11)
#     ax_score.set_title(f'{label} — marker score', fontsize=12)
#     ax_score.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
#     ax_score.grid(axis='y', alpha=0.3)

#     # --- log2fc panel ---
#     ax_log2fc.axhline(0, color='black', linestyle='--', linewidth=1.2,
#                       alpha=0.6, label='Baseline (preflight mean)')
#     ax_log2fc.axhspan(0,  ax_log2fc.get_ylim()[1] if ax_log2fc.get_ylim()[1] > 0 else 2,
#                       alpha=0.04, color='red',  label='Upregulated')
#     ax_log2fc.axhspan(ax_log2fc.get_ylim()[0] if ax_log2fc.get_ylim()[0] < 0 else -2,
#                       0, alpha=0.04, color='blue', label='Downregulated')
#     ax_log2fc.set_ylabel('Mean log₂FC vs preflight baseline', fontsize=11)
#     ax_log2fc.set_title(f'{label} — log₂FC', fontsize=12)
#     ax_log2fc.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
#     ax_log2fc.set_xlabel('Timepoint', fontsize=11)
#     ax_log2fc.tick_params(axis='x', rotation=45)
#     ax_log2fc.grid(axis='y', alpha=0.3)

#     fig.suptitle(
#         f'{label} — Recovery Score and Log₂FC Over Time\n'
#         f'Post-flight vs Preflight Baseline (OSD-575 EVE Panel)',
#         fontsize=13
#     )
#     fig.tight_layout()
#     fig.subplots_adjust(top=0.88)
#     plt.show()
