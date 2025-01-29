import re
import pandas as pd

# Function to clean peptide sequence
def clean_peptide(peptide):
    return re.sub(r'\[.*?\]', '', peptide)

# Read merged.txt file
with open('merged.txt', 'r') as file:
    lines = file.readlines()

# Dictionary to store peptide counts
peptide_counts = {}

# Process each line in merged.txt
for line in lines:
    if 'mzspec:PXD012308:' in line:
        parts = line.split(':')
        peptide_identification = parts[5].split('/')[0]
        cleaned_peptide = clean_peptide(peptide_identification)
        if cleaned_peptide in peptide_counts:
            peptide_counts[cleaned_peptide] += 1
            print(cleaned_peptide, peptide_counts[cleaned_peptide])
        else:
            peptide_counts[cleaned_peptide] = 1

# Read the example_reanalysis._peptide.py file
df = pd.read_csv('example_reanalysis_peptide.tsv', sep='\t')

# Update the PA_psms column
df['PA_psms'] = df['Peptide sequence'].apply(lambda x: peptide_counts.get(clean_peptide(x), 0))

# Write the updated DataFrame back to the file
df.to_csv('example_reanalysis_peptide.tsv', sep='\t', index=False)