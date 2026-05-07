from cytokine_dataframe import *

# Dataset straight from open-source
immune_alamar = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum.immune.AlamarPanel_TRANSFORMED.csv', index_col=0).transpose())

# Clean up dataset
immune_alamar_long = process_cytokine_df(immune_alamar, 'immune_alamar')