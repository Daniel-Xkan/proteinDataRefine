import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import os

# Load the data
chunk_size = 10000
fasta_file = 'uniprotkb_human_proteome_UP000005640_with_isoforms_2024-10-08.fasta'

# Define a function to get the peptide sequence
def get_peptide_sequence(row):
    return row['PeptideAtlas_peptide_demod'] if pd.notna(row['PeptideAtlas_peptide_demod']) else row['opt_global_UnmodPep']

# Define a function to get the peptide charge
def get_peptide_charge(row):
    return row['PeptideAtlas_charge'] if pd.notna(row['PeptideAtlas_charge']) else row['charge']

# Define a function to get the protein ID from the sequence column "tr|D9J307|D9J307_HUMAN"
def get_protein_id_from_msgf(row):
    accession = row['accession']
    if isinstance(accession, str) and '|' in accession:
        parts = accession.split('|')
        if len(parts) > 1:
            return parts[1]
    return None

# Load protein sequences from fasta file
# Pre-compute and store sequences in a dictionary for faster lookups
protein_sequences = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    gene_id = next((part.split('=')[1] for part in record.description.split() if part.startswith('GN=')), 'UNKNOWN')
    protein_sequences[record.id] = {
        'gene_id': gene_id,
        'sequence': str(record.seq).replace('I', 'L')
    }

# Function to count matches in protein sequences
def count_matches(peptide, allow_mutation=False):
    peptide = peptide.replace('I', 'L')
    protein_ids = set()
    gene_ids = set()
    
    for header, data in protein_sequences.items():
        gene_id = data['gene_id']
        sequence = data['sequence']
        protein_id = header.split('|')[1] if '|' in header else header
        
        if allow_mutation:
            # Check for near-matches (SAAP)
            for i in range(len(sequence) - len(peptide) + 1):
                window = sequence[i:i+len(peptide)]
                if sum(1 for a, b in zip(peptide, window) if a != b) <= 1:
                    protein_ids.add(protein_id)
                    gene_ids.add(gene_id)
                    break
        else:
            # Check for exact matches
            if peptide in sequence:
                protein_ids.add(protein_id)
                gene_ids.add(gene_id)
    
    # Return counts and lists
    return (
        len(protein_ids), 
        len(gene_ids), 
        ';'.join(protein_ids) if protein_ids else 'None', 
        ';'.join(gene_ids) if gene_ids else 'UNKNOWN'
    )

# Function to process a peptide
peptideatlas_df = pd.read_csv('PeptideAtlas_peptides.tsv', sep='\t')
def process_peptide(to_parallel_process):
    peptide, peptide_data = to_parallel_process
    print(f"Processing new peptide: {peptide}")
    # Get exact matches and SAAP matches
    num_proteins, num_genes, list_proteins, list_genes = count_matches(peptide)
    num_proteins_saap, num_genes_saap, list_proteins_saap, list_genes_saap = count_matches(peptide, allow_mutation=True)
    
    # Build the output row
    peptide_row = {
        'Peptide sequence': peptide,
        'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
        'Protein identifier': peptide_data.apply(get_protein_id_from_msgf, axis=1).iloc[0],
        'Num_specs_both': 0,
        'Num_specs_MSGF': 0,
        'Num_specs_PA': 0,
        'PA_peptide': 0,
        'PA_psms':0,
        'Num_proteins': num_proteins,
        'List_proteins': list_proteins,
        'Num_genes': num_genes,
        'List_genes': list_genes,
        'Num_proteins_saap': num_proteins_saap,
        'List_proteins_saap': list_proteins_saap,
        'Num_genes_saap': num_genes_saap,
        'List_genes_saap': list_genes_saap
    }
    #print(peptide_row)
    return peptide_row
    # Function to process an existing peptide
###########################################################################
# Load the peptide and spectrum data
#peptide_data_path = "example_reanalysis_peptide.tsv"
def main():
    spectrum_data_path = "example_reanalysis_spectrum.tsv"

    # Load spectrum data
    spectrum_df = pd.read_csv(spectrum_data_path, sep="\t")

    # Create a set of unique peptide IDs with their charge and protein identifier as tuples
    unique_peptides = set(
        (get_peptide_sequence(row), get_peptide_charge(row), get_protein_id_from_msgf(row))
        for _, row in spectrum_df.iterrows()
    )

    peptide_sequences_set = set()
    # Check if the file already exists
    if os.path.exists("example_reanalysis_peptide.tsv"):
        # Load the existing file
        peptide_df = pd.read_csv("example_reanalysis_peptide.tsv", sep="\t")
    else:
        # Create a new DataFrame with the required headers
        headers = [
            'Peptide sequence', 'Peptide charge', 'Protein identifier', 'Num_specs_both',
            'Num_specs_MSGF', 'Num_specs_PA', 'PA_peptide', 'PA_psms', 'Num_proteins',
            'List_proteins', 'Num_genes', 'List_genes', 'Num_proteins_saap',
            'List_proteins_saap', 'Num_genes_saap', 'List_genes_saap'
        ]
        peptide_df = pd.DataFrame(columns=headers)
        # Save the initial peptide dataframe
        peptide_df.to_csv("example_reanalysis_peptide.tsv", sep="\t", index=False)

    # Add all peptide sequences from the peptide dataframe to a set
    peptide_sequences_set.update(peptide_df['Peptide sequence'])

    # Initialize a counter for unique peptides added
    unique_peptides_added_count = 0
    print(f"we have {len(unique_peptides)} new unique peptides.")

    with Pool(cpu_count()) as pool:
        # Use unordered imap to process peptides in parallel
        results = pool.imap_unordered(
            process_peptide,
            ((peptide, spectrum_df[spectrum_df['sequence'] == peptide]) for peptide, _, _ in unique_peptides if peptide not in peptide_sequences_set)
        )
        
        for new_row in results:
            peptide_sequences_set.add(new_row['Peptide sequence'])
            peptide_df = pd.concat([peptide_df, pd.DataFrame([new_row])], ignore_index=True)

    # Save the updated dataframe
    peptide_df.to_csv("example_reanalysis_peptide.tsv", sep="\t", index=False)

    # Function to update the specified columns
    def update_peptide_data(peptide_df, spectrum_df):
        # Iterate through each peptide in the peptide dataframe
        for index, row in peptide_df.iterrows():
            peptide_sequence = row['Peptide sequence']
            
            # Filter spectrum data for the current peptide
            spectrum_subset = spectrum_df[spectrum_df['sequence'] == peptide_sequence]
            
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

if __name__ == "__main__":
    main()
