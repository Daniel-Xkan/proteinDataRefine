import pandas as pd
import re

# Load the peptide and spectrum data
peptide_data_path = "example_peptide_reanalysis.tsv"
spectrum_data_path = "example_reanalysis_spectrum.tsv"
peptideatlas_df = pd.read_csv('PeptideAtlas_peptides.tsv', sep="\t")

peptide_df = pd.read_csv(peptide_data_path, sep="\t")
spectrum_df = pd.read_csv(spectrum_data_path, sep="\t")


def get_protein_id_from_pa(row):
    peptide_sequence = row['Peptide sequence']
    matching_rows = peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide_sequence]
    if not matching_rows.empty:
        return matching_rows.iloc[0, 0]  # Return the protein ID from column 0
    return None
output_df = pd.DataFrame()

# Function to update the specified columns
def update_peptide_data(peptide_df, spectrum_df):
    # Iterate through each peptide in the peptide dataframe
    pa_observations_df = pd.read_csv('PA_observations.csv')
    for index, row in peptide_df.iterrows():
        # Print progress update
        if index % 100 == 0:
            print(f"Processing peptide {index+1}/{len(peptide_df)} ({(index+1)/len(peptide_df)*100:.1f}%)")
        peptide_sequence = row['Peptide sequence']
        
        # Filter spectrum data for the current peptide
        spectrum_subset = spectrum_df[(spectrum_df['sequence'] == peptide_sequence) | (spectrum_df['PeptideAtlas_peptide_demod'] == peptide_sequence)]
        
        # Update the columns
        num_specs_both = len(spectrum_subset[(pd.notna(spectrum_subset['PeptideAtlas_USI'])) & (pd.notna(spectrum_subset['sequence']))])
        num_specs_msgf = len(spectrum_subset[(pd.isna(spectrum_subset['PeptideAtlas_USI'])) & (pd.notna(spectrum_subset['sequence']))])
        num_specs_pa = len(spectrum_subset[(pd.notna(spectrum_subset['PeptideAtlas_USI'])) & (pd.isna(spectrum_subset['sequence']))])
        pa_peptide = 0

        if peptide_sequence in pa_observations_df['Sequence'].values:
            pa_peptide = 1
        # Count PSMs from all_usi.txt
        pa_psms = 0
        # try:
        #     with open('all_usii.txt', 'r') as usi_file:
        #         for line in usi_file:
        #             usi = line.strip()
        #             # Extract and demodify the peptide sequence from the USI
        #             if ":" not in usi or "/" not in usi:
        #                 continue
                    
        #             # Extract sequence part between last colon and slash
        #             parts = usi.split(':')
        #             if len(parts) < 2:
        #                 continue
                    
        #             seq_part = parts[-1].split('/')[0]
                    
        #             # Remove modifications (text in square brackets)
        #             demod_seq = re.sub(r'\[.*?\]', '', seq_part)
                    
        #             # If peptide matches the current peptide sequence
        #             if demod_seq == peptide_sequence:
        #                 pa_psms += 1
        # except FileNotFoundError:
        # pa_psms = len(spectrum_subset[spectrum_subset['sequence'] == peptide_sequence])
        
        # Update the peptide dataframe
        
        protein_identifier = row['Protein identifier']
        if isinstance(protein_identifier, str):
            peptide_df.at[index, 'Protein identifier MSGF'] = protein_identifier.split('-')[0]
        else:
            peptide_df.at[index, 'Protein identifier MSGF'] = None
        peptide_df.at[index, 'Protein identifier PA'] = ";".join(peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide_sequence].iloc[:, 0].dropna().unique())
        peptide_df.at[index, 'Num_specs_both'] = num_specs_both
        peptide_df.at[index, 'Num_specs_MSGF'] = num_specs_msgf
        peptide_df.at[index, 'Num_specs_PA'] = num_specs_pa
        peptide_df.at[index, 'PA_peptide'] = pa_peptide
        peptide_df.at[index, 'PA_psms'] = pa_psms

    return peptide_df

# Update the peptide data
updated_peptide_df = update_peptide_data(peptide_df, spectrum_df)

# Save the updated dataframe
updated_peptide_df.to_csv("example_peptide_reanalysis_updated.tsv", sep="\t", index=False)