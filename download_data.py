"""
download_data.py — Pre-download all NASA OSDR datasets locally
Run once before run_all_scores.py to avoid network dependency during demo.
"""

import os
import pandas as pd

os.makedirs('data', exist_ok=True)

print("Downloading CBC...")
pd.read_csv(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager'
    '&file=LSDS-7_Complete_Blood_Count_CBC.upload_SUBMITTED.csv',
    index_col=0
).to_csv('data/cbc.csv')
print("  Saved: data/cbc.csv")

print("Downloading m6A...")
xl = pd.read_excel(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager'
    '&file=GLDS-561_directm6Aseq_Direct_RNA_seq_m6A_Processed_Data.xlsx',
    sheet_name='I4-FP1',
    skiprows=[0, 1, 2, 3, 4, 5, 6, 7],
    header=[0, 1],
    index_col=0
)
xl.columns = [f"{str(a)}_{str(b)}".strip("_") for a, b in xl.columns]
xl.to_csv('data/m6a_FP1.csv')
print("  Saved: data/m6a_FP1.csv")

print("Downloading cytokine EVE panel...")
pd.read_csv(
    'https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager'
    '&file=LSDS-8_Multiplex_serum_immune_EvePanel_TRANSFORMED.csv',
    index_col=0
).to_csv('data/cytokine_eve.csv')
print("  Saved: data/cytokine_eve.csv")

print("\nDone. All datasets saved to data/")
