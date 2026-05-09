import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# =============================================================================
# CBC ANALYSIS — OSD-569
#
# MD/PhD advice: either use raw CBC values or don't use them.
# We use raw values with clinical reference range bands so the data is
# interpretable in context — a WBC of 11.2 is meaningful; a log2FC of 0.3 is not.
#
# CBC is used here as corroborating evidence for the RNAseq anemia story,
# NOT as a primary scoring domain. The RNAseq erythroid domain carries the
# scientific argument; CBC provides the clinical validation.
#
# NLR (Neutrophil-to-Lymphocyte Ratio) is included as a derived inflammatory
# marker — it's a validated clinical stress index and does not require
# reference ranges to interpret (higher = more stress).
# =============================================================================

# Clinical reference ranges (adult males, standard lab values)
# These are approximate and vary by lab — flag, don't diagnose
REFERENCE_RANGES = {
    'WHITE BLOOD CELL COUNT':  (4.5,  11.0,  '10³/µL'),
    'NEUTROPHILS':             (1.8,  7.7,   '10³/µL'),
    'LYMPHOCYTES':             (1.0,  4.8,   '10³/µL'),
    'RED BLOOD CELL COUNT':    (4.7,  6.1,   '10⁶/µL'),
    'HEMOGLOBIN':              (13.5, 17.5,  'g/dL'),
    'PLATELET COUNT':          (150,  400,   '10³/µL'),
}

TIMEPOINT_ORDER = ['L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82', 'R+194']
PREFLIGHT       = ['L-92', 'L-44', 'L-3']
POSTFLIGHT      = ['R+1', 'R+45', 'R+82', 'R+194']
TIMEPOINT_DAYS  = {
    'L-92': -92, 'L-44': -44, 'L-3': -3,
    'R+1': 1, 'R+45': 45, 'R+82': 82, 'R+194': 194
}

metric_col   = 'ANALYTE'
subject_col  = 'SUBJECT_ID'
timepoint_col= 'TEST_DATE'
value_col    = 'VALUE'


# =============================================================================
# LOAD
# =============================================================================

cbc = pd.read_csv(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download'
    '?source=datamanager&file=LSDS-7_Complete_Blood_Count_CBC.upload_SUBMITTED.csv',
    index_col=0
)

cbc_wide = cbc.pivot_table(
    index=[subject_col, timepoint_col],
    columns=metric_col,
    values=value_col
).reset_index()

# Drop female RBC count — cohort is male
cbc_wide = cbc_wide.drop(columns=['RED BLOOD CELL COUNT (FEMALE)'], errors='ignore')

subjects = sorted(cbc_wide[subject_col].unique())


# =============================================================================
# PLOT 1: RAW VALUES WITH REFERENCE RANGES
# Shows actual clinical values, not fold-changes.
# Reference range shaded in green; out-of-range regions in pink.
# =============================================================================

def plot_cbc_raw(output_path='fig_cbc_raw_values.png'):
    metrics = [m for m in REFERENCE_RANGES if m in cbc_wide.columns]
    n = len(metrics)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    colors = plt.cm.tab10.colors

    for i, metric in enumerate(metrics):
        ax = axes[i]
        low, high, unit = REFERENCE_RANGES[metric]

        # Reference range band
        ax.axhspan(low, high, color='#e8f5e9', alpha=0.6, label='Normal range', zorder=0)
        ax.axhline(low,  color='#81c784', linewidth=0.8, linestyle='--', zorder=1)
        ax.axhline(high, color='#81c784', linewidth=0.8, linestyle='--', zorder=1)

        # Launch line
        ax.axvline(3.5, color='#e74c3c', linewidth=1.2, linestyle='--',
                   alpha=0.5, label='Landing', zorder=2)

        for j, subject in enumerate(subjects):
            sub = cbc_wide[cbc_wide[subject_col] == subject].copy()
            sub = sub[sub[timepoint_col].isin(TIMEPOINT_ORDER)]
            sub['tp_order'] = sub[timepoint_col].map({t: k for k, t in enumerate(TIMEPOINT_ORDER)})
            sub = sub.sort_values('tp_order')

            ax.plot(sub['tp_order'], sub[metric],
                    marker='o', linewidth=1.8, markersize=5,
                    color=colors[j % len(colors)],
                    label=subject, alpha=0.85, zorder=3)

        ax.set_xticks(range(len(TIMEPOINT_ORDER)))
        ax.set_xticklabels(TIMEPOINT_ORDER, rotation=45, ha='right', fontsize=7)
        ax.set_title(metric.title(), fontsize=9, fontweight='bold')
        ax.set_ylabel(unit, fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Legend on last axis
    handles = [mpatches.Patch(color=colors[j], label=s) for j, s in enumerate(subjects)]
    handles.append(mpatches.Patch(color='#e8f5e9', label='Normal range'))
    axes[n - 1].legend(handles=handles, fontsize=8, loc='best')

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.suptitle(
        'CBC Raw Values with Clinical Reference Ranges\n'
        'OSD-569 — Pre-flight through R+194',
        fontsize=13, fontweight='bold'
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f'[✓] Saved: {output_path}')
    plt.close(fig)


# =============================================================================
# PLOT 2: NEUTROPHIL-TO-LYMPHOCYTE RATIO
# NLR > 3 is considered elevated; > 5 is markedly elevated.
# Useful as a corroborating inflammation marker alongside cytokine data.
# =============================================================================

def plot_nlr(output_path='fig_cbc_nlr.png'):
    cbc_nlr = cbc_wide.copy()
    cbc_nlr['NLR'] = cbc_nlr['NEUTROPHILS'] / cbc_nlr['LYMPHOCYTES']
    cbc_nlr['days'] = cbc_nlr[timepoint_col].map(TIMEPOINT_DAYS)
    cbc_nlr = cbc_nlr[cbc_nlr[timepoint_col].isin(TIMEPOINT_ORDER)]

    colors = plt.cm.tab10.colors
    fig, ax = plt.subplots(figsize=(11, 5))

    # Reference bands
    ax.axhspan(0, 3,  color='#e8f5e9', alpha=0.4, label='Normal NLR (< 3)')
    ax.axhspan(3, 5,  color='#fff9c4', alpha=0.4, label='Elevated NLR (3–5)')
    ax.axhspan(5, 15, color='#ffebee', alpha=0.4, label='Markedly elevated (> 5)')

    ax.axvline(0, color='#e74c3c', linestyle='--', linewidth=1.2, alpha=0.6, label='Landing')

    for j, subject in enumerate(subjects):
        sub = cbc_nlr[cbc_nlr[subject_col] == subject].sort_values('days')
        ax.plot(sub['days'], sub['NLR'],
                marker='o', linewidth=2, markersize=6,
                color=colors[j % len(colors)], label=subject)

    ax.set_xticks(list(TIMEPOINT_DAYS.values()))
    ax.set_xticklabels(list(TIMEPOINT_DAYS.keys()), rotation=45, ha='right')
    ax.set_xlabel('Timepoint')
    ax.set_ylabel('Neutrophil-to-Lymphocyte Ratio')
    ax.set_title(
        'Neutrophil-to-Lymphocyte Ratio (NLR) Over Time\n'
        'OSD-569 — Inflammatory stress index',
        fontsize=12, fontweight='bold'
    )
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight', dpi=150)
    print(f'[✓] Saved: {output_path}')
    plt.close(fig)


# =============================================================================
# ANEMIA CORROBORATION TABLE
# Summary of hemoglobin and RBC count vs reference range at each postflight
# timepoint — feeds directly into the anemia narrative.
# =============================================================================

def anemia_summary():
    """
    Print a table showing whether each astronaut's hemoglobin and RBC count
    fall below the lower bound of the clinical reference range at each
    postflight timepoint. Used as supporting evidence for the RNAseq
    erythroid domain score.
    """
    anemia_markers = {
        'HEMOGLOBIN':         REFERENCE_RANGES['HEMOGLOBIN'][0],
        'RED BLOOD CELL COUNT': REFERENCE_RANGES['RED BLOOD CELL COUNT'][0],
    }
    post = cbc_wide[cbc_wide[timepoint_col].isin(POSTFLIGHT)]

    print("\n=== ANEMIA CORROBORATION (below lower reference bound = flagged) ===")
    for metric, lower_bound in anemia_markers.items():
        if metric not in post.columns:
            continue
        print(f"\n{metric}  (lower bound: {lower_bound})")
        pivot = post.pivot_table(index=subject_col, columns=timepoint_col, values=metric)
        # Reorder columns
        ordered = [t for t in POSTFLIGHT if t in pivot.columns]
        pivot = pivot[ordered]
        for col in ordered:
            pivot[col] = pivot[col].apply(
                lambda v: f'{v:.2f} ⚠️' if pd.notna(v) and v < lower_bound else f'{v:.2f}'
            )
        print(pivot.to_string())


# =============================================================================
# RUN
# =============================================================================

if __name__ == '__main__':
    plot_cbc_raw()
    plot_nlr()
    anemia_summary()

# # work specifically on whole blood rna seq
# # function in one file, execution in another, plot in another

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns

# cbc = pd.read_csv(
#     'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager&file=LSDS-7_Complete_Blood_Count_CBC.upload_SUBMITTED.csv',
#     index_col=0
# )

# #print(cbc.head())
# #print(cbc.shape)
# #print(cbc.columns.tolist())

# # Check unique values in each column to identify structure
# # for col in cbc.columns:
# #     print(f"{col}: {cbc[col].nunique()} unique values")
# #     print(cbc[col].unique()[:5])
# #     print()

# metric_col = 'ANALYTE'
# subject_col = 'SUBJECT_ID'
# timepoint_col = 'TEST_DATE'
# value_col = 'VALUE'

# cbc_wide = cbc.pivot_table(
#     index=[subject_col, timepoint_col],
#     columns=metric_col,
#     values=value_col
# ).reset_index()

# #print(cbc_wide.head())

# # Timepoint order for plotting
# timepoint_order = ['L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82', 'R+194']

# # Define which timepoints are pre-flight baseline
# preflight = ['L-92', 'L-44', 'L-3']
# postflight = ['R+1', 'R+45', 'R+82', 'R+194']


# # Get list of CBC metrics (exclude ID columns)
# cbc_metrics = [col for col in cbc_wide.columns 
#                if col not in [subject_col, timepoint_col]]

# subjects = cbc_wide[subject_col].unique()

# log2fc_records = []

# for subject in subjects:
#     subject_data = cbc_wide[cbc_wide[subject_col] == subject]
    
#     # Calculate pre-flight mean for each metric
#     pre = subject_data[subject_data[timepoint_col].isin(preflight)]
#     pre_means = pre[cbc_metrics].mean()
    
#     # Calculate log2FC for each post-flight timepoint
#     post = subject_data[subject_data[timepoint_col].isin(postflight)]
    
#     for _, row in post.iterrows():
#         record = {
#             subject_col: subject,
#             timepoint_col: row[timepoint_col]
#         }
#         for metric in cbc_metrics:
#             if pre_means[metric] > 0 and row[metric] > 0:
#                 record[metric] = np.log2(row[metric] / pre_means[metric])
#             else:
#                 record[metric] = np.nan
#         log2fc_records.append(record)

# log2fc_df = pd.DataFrame(log2fc_records)
# #print("\nLog2FC DataFrame:")
# #print(log2fc_df.head())

# #drop female RBC count for  male astronauts (it shows as NaN)
# log2fc_df = log2fc_df.drop(columns=['RED BLOOD CELL COUNT (FEMALE)'], errors='ignore')
# #print(log2fc_df.columns.tolist())


# metrics_to_plot = [
#     'WHITE BLOOD CELL COUNT',
#     'NEUTROPHILS',
#     'LYMPHOCYTES',
#     'RED BLOOD CELL COUNT',
#     'HEMOGLOBIN',
#     'PLATELET COUNT'
# ]

# fig, axes = plt.subplots(2, 3, figsize=(15, 8))
# axes = axes.flatten()

# for i, metric in enumerate(metrics_to_plot):
#     ax = axes[i]
    
#     for subject in subjects:
#         subject_data = log2fc_df[log2fc_df[subject_col] == subject]
#         subject_data = subject_data.sort_values(timepoint_col, 
#                            key=lambda x: x.map({t: i for i, t in enumerate(postflight)}))
#         ax.plot(subject_data[timepoint_col], subject_data[metric], 
#                 marker='o', label=subject, alpha=0.7)
    
#     ax.axhline(y=0, color='black', linestyle='--', linewidth=1, label='Baseline')
#     ax.legend(fontsize=7)
#     ax.set_title(metric)
#     ax.set_xlabel('Timepoint')
#     ax.set_ylabel('log2FC from baseline')
#     ax.tick_params(axis='x', rotation=45)

# plt.suptitle('CBC Recovery Trajectories (log2FC vs Pre-flight Baseline)', 
#              fontsize=14)  # removed y=1.02
# plt.tight_layout()
# plt.savefig('cbc_trajectories.png', dpi=150, bbox_inches='tight')
# plt.show()


# # step 7 - neutrophil:lymphocyte ratio over time
# timepoint_days = {
#     'L-92': -92, 'L-44': -44, 'L-3': -3,
#     'R+1': 1, 'R+45': 45, 'R+82': 82, 'R+194': 194
# }

# cbc_wide['NLR'] = cbc_wide['NEUTROPHILS'] / cbc_wide['LYMPHOCYTES']

# #fig, ax = plt.subplots(figsize=(8, 5))
# fig, ax = plt.subplots(figsize=(11, 5))


# for subject in subjects:
#     subject_data = cbc_wide[cbc_wide[subject_col] == subject].copy()
#     subject_data['days'] = subject_data[timepoint_col].map(timepoint_days)
#     subject_data = subject_data.sort_values('days')
#     ax.plot(subject_data['days'], subject_data['NLR'],
#             marker='o', label=subject)

# ax.axvline(x=0, color='red', linestyle='--', alpha=0.5, label='Landing')
# ax.set_xticks(list(timepoint_days.values()))
# #ax.set_xticklabels(list(timepoint_days.keys()), rotation=45)
# ax.set_xticklabels(list(timepoint_days.keys()), rotation=60, ha='right')
# ax.set_title('Neutrophil-to-Lymphocyte Ratio Over Time')
# ax.set_xlabel('Timepoint')
# ax.set_ylabel('NLR')
# ax.legend()
# plt.tight_layout()
# plt.savefig('nlr_trajectory.png', dpi=150, bbox_inches='tight')
# plt.show()
