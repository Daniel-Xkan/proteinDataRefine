import sys
import subprocess
import requests
from bs4 import BeautifulSoup


# Function to create the correct PeptideAtlas URL for a given peptide
def create_peptide_url(peptide_id):
    base_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide"
    url = f"{base_url}?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis={peptide_id}&action=QUERY"
    print(f"Created PeptideAtlas URL: {url}")
    return url

# Function to capture the USI from the Spectrum page
def get_usi_from_spectrum(spectrum_url):
    try:
        print(f"Fetching Spectrum URL: {spectrum_url}")
        response = requests.get(spectrum_url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Locate the USI form and extract the value
        usi_input = soup.find('input', {'name': 'USI'})
        if usi_input:
            usi = usi_input['value']
            print(f"USI found: {usi}")
            if usi.startswith('mzspec:PXD') or usi.startswith('mzspec:MSV'):
                return usi
        else:
            print("No USI found in this spectrum.")
    except requests.exceptions.RequestException as e:
        print(f"Error fetching spectrum URL: {spectrum_url}, Error: {e}")
    
    return None

# Function to extract USIs from the individual spectra table
def process_peptide_spectra(protein_id, peptide_id):
    print(f"Processing peptide: {peptide_id} for protein: {protein_id}")
    url = create_peptide_url(peptide_id)
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Parse individual spectra table
        spectra_table = soup.find('table', {'id': 'individual_spectra'})
        if not spectra_table:
            print(f"No spectra table found for peptide: {peptide_id}")
            return
        
        for row in spectra_table.find_all('tr')[1:]:  # Skip the header row
            spectrum_link = row.find('a', href=True)
            if spectrum_link:
                spectrum_url = "https://db.systemsbiology.net" + spectrum_link['href']
                usi = get_usi_from_spectrum(spectrum_url)
                if usi:
                    peptide_identifier = usi.split(':')[-1].split('/')[0]
                    if '[' in peptide_identifier or ']' in peptide_identifier:
                        print(f"Modified peptide identified: {peptide_identifier}")
                        save_usi(protein_id, usi, modified=True)
                    else:
                        print(f"Unmodified peptide identified: {peptide_identifier}")
                        save_usi(protein_id, usi, modified=False)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching peptide URL: {url}, Error: {e}")

# Function to save the USI to the corresponding file
def save_usi(protein_id, usi, modified=False):
    file_name = f"{protein_id}_{'modified' if modified else 'unmodified'}_peptides_usi.txt"
    print(f"Saving USI to file: {file_name}")
    
    with open(file_name, 'a') as file:
        file.write(f"{usi}\n")

# Function to split the peptide file into smaller parts
def split_file(file_path, num_splits=50):
    print(f"Splitting file: {file_path} into {num_splits} parts")
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Skip the header line
    header = lines[0]
    data_lines = lines[1:]
    
    # Calculate the number of lines per split
    split_size = len(data_lines) // num_splits
    
    for i in range(num_splits):
        split_file_name = f"{file_path.split('.txt')[0]}_part_{i + 1}.txt"
        with open(split_file_name, 'w') as split_file:
            split_file.write(header)  # Write the header
            # Write the corresponding data lines
            start_index = i * split_size
            end_index = start_index + split_size if i < num_splits - 1 else len(data_lines)
            split_file.writelines(data_lines[start_index:end_index])
        
        print(f"Created split file: {split_file_name}")

# Function to process the peptide file
def process_peptide_file(file_path):
    print(f"Processing peptide file: {file_path}")
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            protein_id, peptide_id, sequence = line.strip().split('\t')
            print(f"Processing line: {line.strip()}")
            process_peptide_spectra(protein_id, peptide_id)

# Main function to handle splitting and parallel processing
def main(file_path):
    num_splits = 50  # Define how many splits you want
    split_file(file_path, num_splits)  # Split the file into parts

    # Create a list of the split files
    files = [f"{file_path.split('.txt')[0]}_part_{i}.txt" for i in range(1, num_splits + 1)]

    # Launch subprocesses for each split file
    processes = []
    for file in files:
        proc = subprocess.Popen(['python', 'script.py', file])  # Adjust 'script.py' if needed
        processes.append(proc)

    # Wait for all subprocesses to finish
    for proc in processes:
        proc.wait()

    print("All processes finished.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_peptide_file>")
        sys.exit(1)

    peptide_file = sys.argv[1]
    main(peptide_file)

    print("Processing complete and USIs stored.")
