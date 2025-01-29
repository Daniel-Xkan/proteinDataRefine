import pandas as pd

# Load the datasets
reanalysis_df = pd.read_csv('MSGF-PLUS-AMBIGUITY-81a33a88-group_by_spectrum-main.tsv', sep='\t')
proteins_to_skip_file = 'MassIVE-KB_HPP_proteins.tsv'
with open(proteins_to_skip_file, 'r') as file:
    massive_kb_proteins = set(line.strip() for line in file)

# Process the 'opt_global_TopCanonicalProtein' column to extract the protein ID
reanalysis_df['opt_global_TopCanonicalProtein'] = reanalysis_df['opt_global_TopCanonicalProtein'].apply(
    lambda x: x.split('|')[1] if isinstance(x, str) and '|' in x else x
)

# Filter the reanalysis dataframe
filtered_reanalysis_df = reanalysis_df[~reanalysis_df['opt_global_TopCanonicalProtein'].isin(massive_kb_proteins)]

# Save the filtered dataframe to a new file
filtered_reanalysis_df.to_csv('filtered.tsv', sep='\t', index=False)