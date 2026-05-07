# work on this later, after we get at least 1 per domain. 

# 1. IMPORTS
import pandas as pd
import numpy as np

# 2. LOAD DATA
snrnaseq = pd.read_excel('https://osdr.nasa.gov/geode-py/ws/studies/OSD-570/download?source=datamanager&file=GLDS-562_snRNA-Seq_PBMC_Gene_Expression_snRNA-seq_Processed_Data.xlsx', skiprows=[0,1,2,3,4,5,6], index_col=0)
snatacseq = pd.read_excel('https://osdr.nasa.gov/geode-py/ws/studies/OSD-570/download?source=datamanager&file=GLDS-562_snATAC-Seq_PBMC_Chromatin_Accessibility_snATAC-seq_Processed_Data.xlsx', skiprows=[0,1,2,3,4,5,6], index_col=0)
vdj = pd.read_excel('https://osdr.nasa.gov/geode-py/ws/studies/OSD-570/download?source=datamanager&file=GLDS-562_scRNA-Seq_VDJ_Results.xlsx', skiprows=[0,1,2], index_col=0)



# 3. snRNA-seq ANALYSIS
# (recovery dynamics R+45 vs R+1)

# 4. snATAC-seq ANALYSIS  
# (acute flight effect R+1 vs pre-flight)

# 5. VDJ DIVERSITY SCORING
# (per-astronaut clonal diversity → feeds into domain score)

# 6. EXPORT RESULTS
# save scores to CSV for dashboard to pick up
#results.to_csv('../scores/immune_pbmc_scores.csv'

