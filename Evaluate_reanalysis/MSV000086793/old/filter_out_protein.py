import pandas as pd

# Load the datasets in chunks
chunksize = 10000
proteins_to_skip_file = 'MassIVE-KB_HPP_proteins.tsv'
with open(proteins_to_skip_file, 'r') as file:
    massive_kb_proteins = set(line.strip() for line in file)

filtered_chunks = []
for chunk in pd.read_csv('ambiguity_merged.tsv', sep='\t', low_memory=False, error_bad_lines=False, chunksize=chunksize):
    # Process the 'opt_global_TopCanonicalProtein' column to extract the protein ID
    chunk['opt_global_TopCanonicalProtein'] = chunk['opt_global_TopCanonicalProtein'].apply(
        lambda x: x.split('|')[1] if isinstance(x, str) and '|' in x else x
    )
    # Filter the chunk
    filtered_chunk = chunk[~chunk['opt_global_TopCanonicalProtein'].isin(massive_kb_proteins)]
    filtered_chunks.append(filtered_chunk)

# Concatenate all filtered chunks
filtered_reanalysis_df = pd.concat(filtered_chunks)

# Save the filtered dataframe to a new file
filtered_reanalysis_df.to_csv('filtered.tsv', sep='\t', index=False)