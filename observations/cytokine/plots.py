from immune_eve import *
import matplotlib.pyplot as plt

for composite in to_plot:
    # Sort timepoints numerically
    composite['timepoint_num'] = composite['timepoint'].str.extract(r'R\+(\d+)').astype(int)
    composite = composite.sort_values(['astronaut', 'timepoint_num'])

    timepoint_order = composite.drop_duplicates('timepoint').sort_values('timepoint_num')['timepoint'].tolist()

    # Plot
    fig, ax = plt.subplots(figsize=(11, 6))

    for astronaut, data in composite.groupby('astronaut'):
        data = data.set_index('timepoint').reindex(timepoint_order)
        ax.plot(timepoint_order, data['marker_score'], 
                marker='o', linewidth=2, markersize=6, label=astronaut)

    # Reference line at 0 = baseline
    ax.axhline(0, color='black', linestyle='--', linewidth=1.2, alpha=0.6, label='Baseline (preflight mean)')

    # Shading to make positive/negative intuitive
    ax.axhspan(0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 2, 
            alpha=0.04, color='red', label='Upregulated')
    ax.axhspan(ax.get_ylim()[0] if ax.get_ylim()[0] < 0 else -2, 
            0, alpha=0.04, color='blue', label='Downregulated')

    ax.set_xlabel('Timepoint', fontsize=12)
    ax.set_ylabel('Score\n(mean log₂FC)', fontsize=12)
    ax.set_title('Recovery Score Over Time\n(Post-flight vs. Preflight Baseline)', fontsize=13)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.show()
