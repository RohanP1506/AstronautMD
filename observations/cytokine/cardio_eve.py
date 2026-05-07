import numpy as np
import pandas as pd
from cytokine_dataframe import *

# Dataset straight from open-source
cardio_eve = organize_df(pd.read_csv('https://osdr.nasa.gov/geode-py/ws/studies/OSD-575/download?source=datamanager&file=LSDS-8_Multiplex_serum_cardiovascular_EvePanel_TRANSFORMED.csv', index_col=0).transpose())

# Clean up dataset
cardio_eve_long = process_cytokine_df(cardio_eve, 'cardio_eve')
