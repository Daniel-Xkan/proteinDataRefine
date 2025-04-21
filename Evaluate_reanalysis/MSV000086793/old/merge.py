import os
import pandas as pd

# Directory containing the TSV files
directory = 'c:/Users/xuech/Desktop/UCSD/grad/ccmc/proteinDataRefine/Evaluate_reanalysis/MSV000086793/'

# Output file path
output_file = os.path.join(directory, 'ambiguity_merged.tsv')

# Remove the output file if it already exists
if os.path.exists(output_file):
    os.remove(output_file)

# Iterate over files in the directory
for filename in os.listdir(directory):
    if filename.startswith('MSGF-PLUS-AMBIGUITY') and filename.endswith('.tsv'):
        print(f'Processing file: {filename}')
        filepath = os.path.join(directory, filename)
        try:
            df = pd.read_csv(filepath, sep='\t', low_memory=False)
        except pd.errors.EmptyDataError:
            print(f'Skipping empty file: {filename}')
            continue
        # Append to the output file
        df.to_csv(output_file, sep='\t', index=False, mode='a', header=not os.path.exists(output_file))