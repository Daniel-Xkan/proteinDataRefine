import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

# Load the data
spectrum_df = pd.read_csv('example_reanalysis_spectrum.tsv', sep='\t')
peptideatlas_df = pd.read_csv('PeptideAtlas_peptides.tsv', sep='\t')
fasta_file = 'uniprotkb_human_proteome_UP000005640_with_isoforms_2024-10-08.fasta'

# Define a function to get the peptide sequence
def get_peptide_sequence(row):
    return row['PeptideAtlas_peptide_demod'] if pd.notna(row['PeptideAtlas_peptide_demod']) else row['opt_global_UnmodPep']

# Define a function to get the peptide charge
def get_peptide_charge(row):
    return row['PeptideAtlas_charge'] if pd.notna(row['PeptideAtlas_charge']) else row['charge']

    # Define a function to get the protein ID from the sequence column "tr|D9J307|D9J307_HUMAN"
def get_protein_id_from_msgf(row):
#accession ############################################    
#   accession = row['accession']


#     if isinstance(accession, str) and '|' in accession:
#         parts = accession.split('|')
#         if len(parts) > 1:
#             return parts[1]
#     return None
#######################################################
#TopCanonicalProtein #################################
    TopCanonicalProtein = row['opt_global_TopCanonicalProtein']
    return TopCanonicalProtein
# Create a new DataFrame for the output

def get_protein_id_from_pa(row):
    peptide_sequence = row['Peptide sequence']
    matching_rows = peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide_sequence]
    if not matching_rows.empty:
        return matching_rows.iloc[0, 0]  # Return the protein ID from column 0
    return None
output_df = pd.DataFrame()

# Get unique peptides
spectrum_df['Peptide sequence'] = spectrum_df.apply(get_peptide_sequence, axis=1)
unique_peptides = spectrum_df['Peptide sequence'].unique()

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
    protein_list = []
    gene_list = []
    
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
################################################################################################


# Function to count matches in protein sequences
# def count_matches(peptide, allow_mutation=False):
#     peptide = peptide.replace('I', 'L')
#     num_proteins = 0
#     num_genes = 0
#     protein_ids = set()
#     gene_ids = set()
#     for header, sequence in protein_sequences.items():
#         if allow_mutation:
#             for i in range(len(sequence) - len(peptide) + 1):
#                 if sum(1 for a, b in zip(peptide, sequence[i:i+len(peptide)]) if a != b) <= 1:
#                     protein_ids.add(header.split()[0])
#                     gene_id = next((part[3:] for part in header.split() if part.startswith('GN=')), 'UNKNOWN')
#                     gene_ids.add(gene_id)
#                     break
#         else:
#             if peptide in sequence:
#                 protein_ids.add(header.split()[0])
#                 gene_id = next((part[3:] for part in header.split() if part.startswith('GN=')), 'UNKNOWN')
#                 gene_ids.add(gene_id)
#     return protein_ids, gene_ids




# Populate the output DataFrame
output_rows = []
# for peptide in unique_peptides:
#     peptide_data = spectrum_df[spectrum_df['Peptide sequence'] == peptide]
#     protein_ids, gene_ids = count_matches(peptide)
#     num_proteins = len(protein_ids)
#     num_genes = len(gene_ids)
    
#     protein_ids_saap, gene_ids_saap = count_matches(peptide, allow_mutation=True)
#     num_proteins_saap = len(protein_ids_saap)
#     num_genes_saap = len(gene_ids_saap)
#     peptide_row = {
#         'Peptide sequence': peptide,
#         'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
#         'Protein identifier': peptide_data.apply(get_protein_id_from_msgf, axis=1).iloc[0],
#         'Num_specs_both': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
#         'Num_specs_MSGF': len(peptide_data[(pd.isna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
#         'Num_specs_PA': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.isna(peptide_data['sequence']))]),
#         'PA_peptide': 1 if not peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].empty else 0,
#         'PA_psms': len(peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide]),
#         'Num_proteins': num_proteins,
#         'List_proteins': ';'.join(protein_ids),
#         'Num_genes': num_genes,
#         'List_genes': ';'.join(gene_ids),
#         'Num_proteins_saap': num_proteins_saap,
#         'Num_genes_saap': num_genes_saap
#     }
#     output_rows.append(peptide_row)

#deepseek version:####################################################
def process_peptide(peptide):
    peptide_data = spectrum_df[spectrum_df['Peptide sequence'] == peptide]
    # Get exact matches and SAAP matches
    num_proteins, num_genes, list_proteins, list_genes = count_matches(peptide)
    num_proteins_saap, num_genes_saap, list_proteins_saap, list_genes_saap = count_matches(peptide, allow_mutation=True)
    
    # Build the output row
    peptide_row = {
        'Peptide sequence': peptide,
        'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
        'Protein identifier MSGF': peptide_data.apply(get_protein_id_from_msgf, axis=1).iloc[0],
        'Protein identifier PA': peptide_data.apply(get_protein_id_from_pa, axis=1).iloc[0],
        'Num_specs_both': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
        'Num_specs_MSGF': len(peptide_data[(pd.isna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
        'Num_specs_PA': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.isna(peptide_data['sequence']))]),
        'PA_peptide': 1 if not peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].empty else 0,
        'PA_psms': len(peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide]),
        'Num_proteins': num_proteins,
        'List_proteins': list_proteins,
        'Num_genes': num_genes,
        'List_genes': list_genes,
        'Num_proteins_saap': num_proteins_saap,
        'List_proteins_saap': list_proteins_saap,
        'Num_genes_saap': num_genes_saap,
        'List_genes_saap': list_genes_saap
    }
    print(peptide_row)
    return peptide_row

# Use multiprocessing to process peptides in parallel
if __name__ == '__main__':
    counter = 0
    total_peptides = len(unique_peptides)
    
    def update_counter(result):
        global counter
        counter += 1
        print(f"Processed {counter}/{total_peptides} peptides")

    with Pool(cpu_count()) as pool:
        output_rows = []
        for peptide_row in pool.imap_unordered(process_peptide, unique_peptides):
            output_rows.append(peptide_row)
            update_counter(peptide_row)

    output_df = pd.DataFrame(output_rows)

    # Save the output DataFrame to a TSV file
    output_df.to_csv('example_reanalysis_peptide.tsv', sep='\t', index=False)