import os
import pandas as pd

# Define the folder path
folder_path = 'c:/Users/xuech/Desktop/UCSD/grad/ccmc/proteinDataRefine/Evaluate_reanalysis/MSV000086793/'

# Initialize an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

# Iterate over all files in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('_protein_reanalysis.tsv'):
        file_path = os.path.join(folder_path, file_name)
        # Read the current file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        if merged_df.empty:
            merged_df = df
        else:
            # Merge the current DataFrame with the merged DataFrame
            merged_df = merged_df.merge(df, on='Protein', suffixes=('', '_dup'), how='outer')
            # Sum the columns with the same name
            for col in df.columns[1:]:
                merged_df[col] = merged_df[col].fillna(0) + merged_df[col + '_dup'].fillna(0)
                merged_df.drop(columns=[col + '_dup'], inplace=True)

# Save the merged DataFrame to a new file
merged_df.to_csv(os.path.join(folder_path, 'merged_protein.tsv'), sep='\t', index=False)