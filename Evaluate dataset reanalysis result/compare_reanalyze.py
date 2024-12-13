import pandas as pd

# Read the PeptideAtlas USIs from merged.txt
usi_file = 'merged.txt'
usi_data = []
with open(usi_file, 'r') as file:
    for line in file:
        parts = line.strip().split(':')
        if len(parts) == 6:
            dataset, spectrum_file, scan, peptide_identification, peptide_charge = parts[1], parts[2], parts[4], parts[5].split('/')[0], parts[5].split('/')[1]
            usi_data.append([line.strip(), dataset, spectrum_file, scan, peptide_identification, peptide_charge])

usi_df = pd.DataFrame(usi_data, columns=['USI', 'Dataset', 'Spectrum_File', 'Scan_Number', 'Peptide_Identification', 'Peptide_Charge'])

# Create a dictionary for quick lookup
usi_dict = {}
for _, row in usi_df.iterrows():
    key = (row['Spectrum_File'], row['Scan_Number'])
    usi_dict[key] = row

# Read the reanalysis results from MSGF-PLUS-AMBIGUITY-81a33a88-group_by_spectrum-main.tsv
reanalysis_file = 'MSGF-PLUS-AMBIGUITY-81a33a88-group_by_spectrum-main.tsv'
reanalysis_df = pd.read_csv(reanalysis_file, sep='\t')

# Initialize columns for the output DataFrame
reanalysis_df.insert(0, 'PeptideAtlas_USI', '')
reanalysis_df.insert(1, 'PeptideAtlas_peptide', '')
reanalysis_df.insert(2, 'PeptideAtlas_charge', '')

# Match spectra from the PeptideAtlas USIs lists to the reanalysis results
for index, row in reanalysis_df.iterrows():
    original_filepath = row['opt_global_OriginalFilepath']
    scan_number = str(row['opt_global_scan'])
    
    # Extract the spectrum file name from the original file path and remove the .mzML extension if present
    spectrum_file = original_filepath.split('/')[-1].replace('.mzML', '')
    
    key = (spectrum_file, scan_number)
    if key in usi_dict:
        usi_row = usi_dict[key]
        reanalysis_df.at[index, 'PeptideAtlas_USI'] = usi_row['USI']
        reanalysis_df.at[index, 'PeptideAtlas_peptide'] = usi_row['Peptide_Identification']
        reanalysis_df.at[index, 'PeptideAtlas_charge'] = usi_row['Peptide_Charge']
        print(f'Matched file: {usi_row["Spectrum_File"]}, Scan number: {usi_row["Scan_Number"]}')
    
    # if index % 100 == 0:
    #     print(f'Processed {index} rows.')

# Save the results to MSV000088387_reanalysis_compare.tsv
output_file = 'MSV000088387_reanalysis_compare.tsv'
reanalysis_df.to_csv(output_file, sep='\t', index=False)