import pandas as pd
from Bio import SeqIO

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
    accession = row['accession']
    if isinstance(accession, str) and '|' in accession:
        parts = accession.split('|')
        if len(parts) > 1:
            return parts[1]
    return None

# Create a new DataFrame for the output
output_df = pd.DataFrame()

# Get unique peptides
spectrum_df['Peptide sequence'] = spectrum_df.apply(get_peptide_sequence, axis=1)
unique_peptides = spectrum_df['Peptide sequence'].unique()

# Load protein sequences from fasta file
protein_sequences = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    gene_id = next((part.split('=')[1] for part in record.description.split() if part.startswith('GN=')), 'UNKNOWN')
    protein_sequences[record.id] = (gene_id, str(record.seq).replace('I', 'L'))

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

#deepseek version:####################################################
def count_matches(peptide, allow_mutation=False):
    peptide = peptide.replace('I', 'L')
    protein_ids = set()
    gene_ids = set()
    protein_list = []
    gene_list = []
    
    for header, gene_and_sequence in protein_sequences.items():
        gene_id, sequence = gene_and_sequence
        # print(protein_sequences.items())
        header_parts = header.split()
        # Extract protein ID (middle part of the first segment, e.g., "B4E2Q0")
        protein_id = header_parts[0].split('|')[1] if '|' in header_parts[0] else header_parts[0]
        # Extract gene ID (from GN= tag)
        # print(header_parts)
        # gene_id = next((part.split('=')[1] for part in header_parts if part.startswith('GN=')), 'UNKNOWN')
        # print(gene_id)

        
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
        # print(protein_ids)
        # print(gene_ids)
    # Return counts and lists
    return (
        len(protein_ids), 
        len(gene_ids), 
        ';'.join(protein_ids) if protein_ids else 'None', 
        ';'.join(gene_ids) if gene_ids else 'UNKNOWN'
    )
################################################################################################

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
for peptide in unique_peptides:
    peptide_data = spectrum_df[spectrum_df['Peptide sequence'] == peptide]
    # Get exact matches and SAAP matches
    num_proteins, num_genes, list_proteins, list_genes = count_matches(peptide)
    num_proteins_saap, num_genes_saap, list_proteins_saap, list_genes_saap = count_matches(peptide, allow_mutation=True)
    
    # Build the output row
    peptide_row = {
        'Peptide sequence': peptide,
        'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
        'Protein identifier': peptide_data.apply(get_protein_id_from_msgf, axis=1).iloc[0],
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
    output_rows.append(peptide_row)
################################################################################################
    print(peptide_row)
    

output_df = pd.DataFrame(output_rows)

# Save the output DataFrame to a TSV file
output_df.to_csv('example_reanalysis_peptide.tsv', sep='\t', index=False)