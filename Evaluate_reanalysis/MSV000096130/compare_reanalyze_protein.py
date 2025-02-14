import pandas as pd

# Load the data
peptide_file = 'example_peptide_reanalysis.tsv'
protein_file = 'PeptideAtlas_proteins_not_in_MassIVE.tsv'

peptide_df = pd.read_csv(peptide_file, sep='\t')
protein_df = pd.read_csv(protein_file, sep='\t')

# Initialize the result dataframe
result_columns = [
    'Protein', 'Num_peptides_both', 'Num_peptides_MSGF', 'Num_peptides_PA',
    'Num_unique_both', 'Num_unique_MSGF', 'Num_unique_PA',
    'Num_specs_both', 'Num_specs_MSGF', 'Num_specs_PA'
]
result_df = pd.DataFrame(columns=result_columns)

# Iterate over each protein
result_list = []
for protein in protein_df['nextprot_accession']:
    # Filter peptides that match the current protein
    matched_peptides = peptide_df[peptide_df['List_proteins'].str.contains(protein, na=False)]
    
    # Calculate the required counts and sums
    num_peptides_both = len(matched_peptides[matched_peptides['Num_specs_both'] > 0])
    num_peptides_MSGF = len(matched_peptides[matched_peptides['Num_specs_MSGF'] > 0])
    num_peptides_PA = len(matched_peptides[matched_peptides['Num_specs_PA'] > 0])
    
    num_unique_both = len(matched_peptides[(matched_peptides['Num_specs_both'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
    num_unique_MSGF = len(matched_peptides[(matched_peptides['Num_specs_MSGF'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
    num_unique_PA = len(matched_peptides[(matched_peptides['Num_specs_PA'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
    
    num_specs_both = matched_peptides['Num_specs_both'].sum()
    num_specs_MSGF = matched_peptides['Num_specs_MSGF'].sum()
    num_specs_PA = matched_peptides['Num_specs_PA'].sum()
    
    # Append the results to the result list
    result_list.append({
        'Protein': protein,
        'Num_peptides_both': num_peptides_both,
        'Num_peptides_MSGF': num_peptides_MSGF,
        'Num_peptides_PA': num_peptides_PA,
        'Num_unique_both': num_unique_both,
        'Num_unique_MSGF': num_unique_MSGF,
        'Num_unique_PA': num_unique_PA,
        'Num_specs_both': num_specs_both,
        'Num_specs_MSGF': num_specs_MSGF,
        'Num_specs_PA': num_specs_PA
    })

# Convert the result list to a dataframe
result_df = pd.concat([pd.DataFrame([result]) for result in result_list], ignore_index=True)

# Save the result to a new file in TSV format
result_df.to_csv('example_protein_reanalysis.tsv', sep='\t', index=False)