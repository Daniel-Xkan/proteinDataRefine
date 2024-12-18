import pandas as pd
import re

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
reanalysis_file = 'filtered.tsv'
reanalysis_df = pd.read_csv(reanalysis_file, sep='\t')

# Initialize columns for the output DataFrame
reanalysis_df.insert(0, 'PeptideAtlas_USI', '')
reanalysis_df.insert(1, 'PeptideAtlas_peptide', '')
reanalysis_df.insert(2, 'PeptideAtlas_peptide_demod', '')
reanalysis_df.insert(3, 'Peptide_match', '')
reanalysis_df.insert(4, 'PeptideAtlas_charge', '')


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
        # Remove all substrings like "[*]" from PeptideAtlas_peptide
        peptide_demod = re.sub(r'\[.*?\]', '', usi_row['Peptide_Identification'])
        reanalysis_df.at[index, 'PeptideAtlas_peptide_demod'] = peptide_demod
        
        # Set Peptide_match to 1 if PeptideAtlas_peptide_demod matches opt_global_UnmodPep, otherwise 0
        reanalysis_df.at[index, 'Peptide_match'] = 1 if peptide_demod == row['opt_global_UnmodPep'] else 0
        reanalysis_df.at[index, 'PeptideAtlas_charge'] = usi_row['Peptide_Charge']
        print(f'Matched file: {usi_row["Spectrum_File"]}, Scan number: {usi_row["Scan_Number"]}')
    
    # if index % 100 == 0:
    #     print(f'Processed {index} rows.')

    # Identify the datasets and spectrum files that have at least one match
    # matched_datasets = set()
matched_spectrum_files = set()
matched_scans = set()

for index, row in reanalysis_df.iterrows():
    if row['PeptideAtlas_USI']:
        usi_parts = row['PeptideAtlas_USI'].split(':')
        if len(usi_parts) == 6:
            dataset = usi_parts[1].replace('.mzML', '')
            spectrum_file = usi_parts[2]
            matched_spectrum_files.add(spectrum_file)
            matched_scans.add((spectrum_file, usi_parts[4]))

# Add empty USIs for unmatched spectra
for key, usi_row in usi_dict.items():
    spectrum_file, scan_number = key
    if usi_row['Dataset'] == 'PXD012308' and spectrum_file in matched_spectrum_files and key not in matched_scans:
        new_row = {
            'PeptideAtlas_USI': usi_row['USI'],
            'PeptideAtlas_peptide': usi_row['Peptide_Identification'],
            'PeptideAtlas_peptide_demod': re.sub(r'\[.*?\]', '', usi_row['Peptide_Identification']),
            'Peptide_match': 0,
            'PeptideAtlas_charge': usi_row['Peptide_Charge'],
            'opt_global_OriginalFilepath': '',
            'opt_global_scan': '',
            'opt_global_UnmodPep': ''
        }
        reanalysis_df = reanalysis_df.append(new_row, ignore_index=True)


# Save the results to MSV000088387_reanalysis_compare.tsv
output_file = 'example_reanalysis_compare.tsv'
reanalysis_df.to_csv(output_file, sep='\t', index=False)