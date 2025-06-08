import os
import pandas as pd

# Directory containing the TSV files
directory = 'c:/Users/xuech/Desktop/UCSD/grad/ccmc/proteinDataRefine/Evaluate_reanalysis/MSV000096130/'

# List to hold dataframes
dataframes = []

# Iterate over files in the directory
for filename in os.listdir(directory):
    if filename.startswith('MSGF-PLUS-AMBIGUITY') and filename.endswith('.tsv'):
        print(f'Processing file: {filename}')
        filepath = os.path.join(directory, filename)
        try:
            df = pd.read_csv(filepath, sep='\t', low_memory=False)
            dataframes.append(df)
        except pd.errors.EmptyDataError:
            print(f'Skipping empty file: {filename}')
        dataframes.append(df)

# Concatenate all dataframes
merged_df = pd.concat(dataframes, ignore_index=True)

# Save the merged dataframe to a TSV file
merged_df.to_csv(os.path.join(directory, 'ambiguity_merged.tsv'), sep='\t', index=False)