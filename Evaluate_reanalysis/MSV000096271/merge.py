import os
import shutil
import pandas as pd

# Directory containing the TSV files
directory = 'c:/Users/xuech/Desktop/UCSD/grad/ccmc/proteinDataRefine/Evaluate_reanalysis/MSV000096271/'

# Output file path
output_file = os.path.join(directory, 'ambiguity_merged.tsv')

# Archive directory
archive_directory = os.path.join(directory, 'ambiguity_archive')

# Create the archive directory if it doesn't exist
if not os.path.exists(archive_directory):
    os.makedirs(archive_directory)

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
        # Move the processed file to the archive directory
        shutil.move(filepath, os.path.join(archive_directory, filename))