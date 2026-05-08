# work specifically on whole blood rna seq
# function in one file, execution in another, plot in another

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cbc = pd.read_csv(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager&file=LSDS-7_Complete_Blood_Count_CBC.upload_SUBMITTED.csv',
    index_col=0
)

#print(cbc.head())
#print(cbc.shape)
#print(cbc.columns.tolist())

# Check unique values in each column to identify structure
# for col in cbc.columns:
#     print(f"{col}: {cbc[col].nunique()} unique values")
#     print(cbc[col].unique()[:5])
#     print()

metric_col = 'ANALYTE'
subject_col = 'SUBJECT_ID'
timepoint_col = 'TEST_DATE'
value_col = 'VALUE'

cbc_wide = cbc.pivot_table(
    index=[subject_col, timepoint_col],
    columns=metric_col,
    values=value_col
).reset_index()

#print(cbc_wide.head())

# Timepoint order for plotting
timepoint_order = ['L-92', 'L-44', 'L-3', 'R+1', 'R+45', 'R+82', 'R+194']

# Define which timepoints are pre-flight baseline
preflight = ['L-92', 'L-44', 'L-3']
postflight = ['R+1', 'R+45', 'R+82', 'R+194']


# Get list of CBC metrics (exclude ID columns)
cbc_metrics = [col for col in cbc_wide.columns 
               if col not in [subject_col, timepoint_col]]

subjects = cbc_wide[subject_col].unique()

log2fc_records = []

for subject in subjects:
    subject_data = cbc_wide[cbc_wide[subject_col] == subject]
    
    # Calculate pre-flight mean for each metric
    pre = subject_data[subject_data[timepoint_col].isin(preflight)]
    pre_means = pre[cbc_metrics].mean()
    
    # Calculate log2FC for each post-flight timepoint
    post = subject_data[subject_data[timepoint_col].isin(postflight)]
    
    for _, row in post.iterrows():
        record = {
            subject_col: subject,
            timepoint_col: row[timepoint_col]
        }
        for metric in cbc_metrics:
            if pre_means[metric] > 0 and row[metric] > 0:
                record[metric] = np.log2(row[metric] / pre_means[metric])
            else:
                record[metric] = np.nan
        log2fc_records.append(record)

log2fc_df = pd.DataFrame(log2fc_records)
#print("\nLog2FC DataFrame:")
#print(log2fc_df.head())

#drop female RBC count for  male astronauts (it shows as NaN)
log2fc_df = log2fc_df.drop(columns=['RED BLOOD CELL COUNT (FEMALE)'], errors='ignore')
#print(log2fc_df.columns.tolist())


metrics_to_plot = [
    'WHITE BLOOD CELL COUNT',
    'NEUTROPHILS',
    'LYMPHOCYTES',
    'RED BLOOD CELL COUNT',
    'HEMOGLOBIN',
    'PLATELET COUNT'
]

fig, axes = plt.subplots(2, 3, figsize=(15, 8))
axes = axes.flatten()

for i, metric in enumerate(metrics_to_plot):
    ax = axes[i]
    
    for subject in subjects:
        subject_data = log2fc_df[log2fc_df[subject_col] == subject]
        subject_data = subject_data.sort_values(timepoint_col, 
                           key=lambda x: x.map({t: i for i, t in enumerate(postflight)}))
        ax.plot(subject_data[timepoint_col], subject_data[metric], 
                marker='o', label=subject, alpha=0.7)
    
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, label='Baseline')
    ax.legend(fontsize=7)
    ax.set_title(metric)
    ax.set_xlabel('Timepoint')
    ax.set_ylabel('log2FC from baseline')
    ax.tick_params(axis='x', rotation=45)

plt.suptitle('CBC Recovery Trajectories (log2FC vs Pre-flight Baseline)', 
             fontsize=14)  # removed y=1.02
plt.tight_layout()
plt.savefig('cbc_trajectories.png', dpi=150, bbox_inches='tight')
plt.show()


# step 7 - neutrophil:lymphocyte ratio over time
timepoint_days = {
    'L-92': -92, 'L-44': -44, 'L-3': -3,
    'R+1': 1, 'R+45': 45, 'R+82': 82, 'R+194': 194
}

cbc_wide['NLR'] = cbc_wide['NEUTROPHILS'] / cbc_wide['LYMPHOCYTES']

#fig, ax = plt.subplots(figsize=(8, 5))
fig, ax = plt.subplots(figsize=(11, 5))


for subject in subjects:
    subject_data = cbc_wide[cbc_wide[subject_col] == subject].copy()
    subject_data['days'] = subject_data[timepoint_col].map(timepoint_days)
    subject_data = subject_data.sort_values('days')
    ax.plot(subject_data['days'], subject_data['NLR'],
            marker='o', label=subject)

ax.axvline(x=0, color='red', linestyle='--', alpha=0.5, label='Landing')
ax.set_xticks(list(timepoint_days.values()))
#ax.set_xticklabels(list(timepoint_days.keys()), rotation=45)
ax.set_xticklabels(list(timepoint_days.keys()), rotation=60, ha='right')
ax.set_title('Neutrophil-to-Lymphocyte Ratio Over Time')
ax.set_xlabel('Timepoint')
ax.set_ylabel('NLR')
ax.legend()
plt.tight_layout()
plt.savefig('nlr_trajectory.png', dpi=150, bbox_inches='tight')
plt.show()

