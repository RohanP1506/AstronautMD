import pandas as pd
rna_blood = pd.read_excel(
    "https://osdr.nasa.gov/geode-py/ws/studies/OSD-569/download?source=datamanager&file=GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx", skiprows=[0,1,2,3,4,5,6,9], header=[0,1], index_col=0)

print(rna_blood)

# work specifically on whole blood rna seq
# function in one file, execution in another, plot in another
