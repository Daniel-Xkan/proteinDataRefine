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

# Create a new DataFrame for the output
output_df = pd.DataFrame()

# Get unique peptides
spectrum_df['Peptide sequence'] = spectrum_df.apply(get_peptide_sequence, axis=1)
unique_peptides = spectrum_df['Peptide sequence'].unique()

# Load protein sequences from fasta file
protein_sequences = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    protein_sequences[record.id] = str(record.seq).replace('I', 'L')

# Function to count matches in protein sequences
def count_matches(peptide, allow_mutation=False):
    peptide = peptide.replace('I', 'L')
    num_proteins = 0
    num_genes = 0
    protein_ids = set()
    gene_ids = set()
    for header, sequence in protein_sequences.items():
        if allow_mutation:
            for i in range(len(sequence) - len(peptide) + 1):
                if sum(1 for a, b in zip(peptide, sequence[i:i+len(peptide)]) if a != b) <= 1:
                    protein_ids.add(header.split()[0])
                    gene_id = next((part[3:] for part in header.split() if part.startswith('GN=')), 'UNKNOWN')
                    gene_ids.add(gene_id)
                    break
        else:
            if peptide in sequence:
                protein_ids.add(header.split()[0])
                gene_id = next((part[3:] for part in header.split() if part.startswith('GN=')), 'UNKNOWN')
                gene_ids.add(gene_id)
    return len(protein_ids), len(gene_ids)

# Populate the output DataFrame
output_rows = []
for peptide in unique_peptides:
    peptide_data = spectrum_df[spectrum_df['Peptide sequence'] == peptide]
    num_proteins, num_genes = count_matches(peptide)
    num_proteins_saap, num_genes_saap = count_matches(peptide, allow_mutation=True)
    peptide_row = {
        'Peptide sequence': peptide,
        'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
        'Protein identifier': peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].iloc[:, 0].values[0] if not peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].empty else None,
        'Num_specs_both': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
        'Num_specs_MSGF': len(peptide_data[(pd.isna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
        'Num_specs_PA': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.isna(peptide_data['sequence']))]),
        'PA_peptide': 1 if not peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].empty else 0,
        'PA_psms': len(peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide]),
        'Num_proteins': num_proteins,
        'Num_genes': num_genes,
        'Num_proteins_saap': num_proteins_saap,
        'Num_genes_saap': num_genes_saap
    }
    output_rows.append(peptide_row)
    print(f'Processed peptide: {peptide}, Num proteins: {num_proteins}, Num proteins with mutation: {num_proteins_saap}, Num genes: {num_genes}, Num genes with mutation: {num_genes_saap}')
    

output_df = pd.DataFrame(output_rows)

# Save the output DataFrame to a TSV file
output_df.to_csv('example_reanalysis_peptide.tsv', sep='\t', index=False)