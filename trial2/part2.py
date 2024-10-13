
import requests
from bs4 import BeautifulSoup


# Function to create the correct PeptideAtlas URL for a given peptide
def create_peptide_url(peptide_id):
    base_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide"
    return f"{base_url}?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis={peptide_id}&action=QUERY"

# Function to capture the USI from the Spectrum page
def get_usi_from_spectrum(spectrum_url):
    response = requests.get(spectrum_url)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Locate the USI form and extract the value
    usi_input = soup.find('input', {'name': 'USI'})
    if usi_input:
        usi = usi_input['value']
        if usi.startswith('mzspec:PXD') or usi.startswith('mzspec:MSV'):
            return usi
    return None

# Function to extract USIs from the individual spectra table
def process_peptide_spectra(protein_id, peptide_id):
    url = create_peptide_url(peptide_id)
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Parse individual spectra table
    spectra_table = soup.find('table', {'id': 'individual_spectra'})
    if not spectra_table:
        return
    
    for row in spectra_table.find_all('tr')[1:]:  # Skip the header row
        spectrum_link = row.find('a', href=True)
        if spectrum_link:
            spectrum_url = "https://db.systemsbiology.net" + spectrum_link['href']
            usi = get_usi_from_spectrum(spectrum_url)
            if usi:
                peptide_identifier = usi.split(':')[-1].split('/')[0]
                if '[' in peptide_identifier or ']' in peptide_identifier:
                    save_usi(protein_id, usi, modified=True)
                else:
                    save_usi(protein_id, usi, modified=False)

# Function to save the USI to the corresponding file
def save_usi(protein_id, usi, modified=False):
    if modified:
        file_name = f"{protein_id}_modified_peptides_usi.txt"
    else:
        file_name = f"{protein_id}_unmodified_peptides_usi.txt"
    
    with open(file_name, 'a') as file:
        file.write(f"{usi}\n")

# Process the PA2024_HPP_peptides.txt file
def process_peptide_file(file_path):
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            protein_id, peptide_id, sequence = line.strip().split('\t')
            process_peptide_spectra(protein_id, peptide_id)

# Example usage
file_path = 'PA2024_HPP_peptides.txt'
process_peptide_file(file_path)


print("Processing complete and USIs stored.")