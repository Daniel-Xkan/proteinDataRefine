import pandas as pd

# Load the peptide and spectrum data
peptide_data_path = "example_reanalysis_peptide.tsv"
spectrum_data_path = "example_reanalysis_spectrum.tsv"

peptide_df = pd.read_csv(peptide_data_path, sep="\t")
spectrum_df = pd.read_csv(spectrum_data_path, sep="\t")

# Function to update the specified columns
def update_peptide_data(peptide_df, spectrum_df):
    # Iterate through each peptide in the peptide dataframe
    for index, row in peptide_df.iterrows():
        peptide_sequence = row['Peptide sequence']
        
        # Filter spectrum data for the current peptide
        spectrum_subset = spectrum_df[(spectrum_df['sequence'] == peptide_sequence) | (spectrum_df['PeptideAtlas_peptide_demod'] == peptide_sequence)]
        
        # Update the columns
        num_specs_both = len(spectrum_subset[(pd.notna(spectrum_subset['PeptideAtlas_USI'])) & (pd.notna(spectrum_subset['sequence']))])
        num_specs_msgf = len(spectrum_subset[(pd.isna(spectrum_subset['PeptideAtlas_USI'])) & (pd.notna(spectrum_subset['sequence']))])
        num_specs_pa = len(spectrum_subset[(pd.notna(spectrum_subset['PeptideAtlas_USI'])) & (pd.isna(spectrum_subset['sequence']))])
        pa_peptide = 1 if not spectrum_subset[spectrum_subset['sequence'] == peptide_sequence].empty else 0
        pa_psms = len(spectrum_subset[spectrum_subset['sequence'] == peptide_sequence])
        
        # Update the peptide dataframe
        peptide_df.at[index, 'Num_specs_both'] = num_specs_both
        peptide_df.at[index, 'Num_specs_MSGF'] = num_specs_msgf
        peptide_df.at[index, 'Num_specs_PA'] = num_specs_pa
        peptide_df.at[index, 'PA_peptide'] = pa_peptide
        peptide_df.at[index, 'PA_psms'] = pa_psms

    return peptide_df

# Update the peptide data
updated_peptide_df = update_peptide_data(peptide_df, spectrum_df)

# Save the updated dataframe
updated_peptide_df.to_csv("example_reanalysis_peptide.tsv", sep="\t", index=False)