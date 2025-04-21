import pandas as pd
import re
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import os

# filter_spectrum_peptide_protein.py

def compare_peptide(spectrum_file, peptideatlas_file, fasta_file, output_file):
    spectrum_df = pd.read_csv(spectrum_file, sep='\t')
    peptideatlas_df = pd.read_csv(peptideatlas_file, sep='\t')

    def get_peptide_sequence(row):
        return row['PeptideAtlas_peptide_demod'] if pd.notna(row['PeptideAtlas_peptide_demod']) else row['opt_global_UnmodPep']

    def get_peptide_charge(row):
        return row['PeptideAtlas_charge'] if pd.notna(row['PeptideAtlas_charge']) else row['charge']

    def get_protein_id_from_msgf(row):
        accession = row['accession']
        if isinstance(accession, str) and '|' in accession:
            parts = accession.split('|')
            if len(parts) > 1:
                return parts[1]
        return None

    spectrum_df['Peptide sequence'] = spectrum_df.apply(get_peptide_sequence, axis=1)
    unique_peptides = spectrum_df['Peptide sequence'].unique()

    protein_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = next((part.split('=')[1] for part in record.description.split() if part.startswith('GN=')), 'UNKNOWN')
        protein_sequences[record.id] = {
            'gene_id': gene_id,
            'sequence': str(record.seq).replace('I', 'L')
        }

    def count_matches(peptide, allow_mutation=False):
        peptide = peptide.replace('I', 'L')
        protein_ids = set()
        gene_ids = set()
        for header, data in protein_sequences.items():
            gene_id = data['gene_id']
            sequence = data['sequence']
            protein_id = header.split('|')[1] if '|' in header else header
            if allow_mutation:
                for i in range(len(sequence) - len(peptide) + 1):
                    window = sequence[i:i+len(peptide)]
                    if sum(1 for a, b in zip(peptide, window) if a != b) <= 1:
                        protein_ids.add(protein_id)
                        gene_ids.add(gene_id)
                        break
            else:
                if peptide in sequence:
                    protein_ids.add(protein_id)
                    gene_ids.add(gene_id)
        return (
            len(protein_ids), 
            len(gene_ids), 
            ';'.join(protein_ids) if protein_ids else 'None', 
            ';'.join(gene_ids) if gene_ids else 'UNKNOWN'
        )

    def process_peptide(peptide):
        peptide_data = spectrum_df[spectrum_df['Peptide sequence'] == peptide]
        num_proteins, num_genes, list_proteins, list_genes = count_matches(peptide)
        num_proteins_saap, num_genes_saap, list_proteins_saap, list_genes_saap = count_matches(peptide, allow_mutation=True)
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
        return peptide_row

    if __name__ == '__main__':
        with Pool(cpu_count()) as pool:
            output_rows = []
            for peptide_row in pool.imap_unordered(process_peptide, unique_peptides):
                output_rows.append(peptide_row)

        output_df = pd.DataFrame(output_rows)
        output_df.to_csv(output_file, sep='\t', index=False)

def compare_protein(peptide_file, protein_file, output_file):
    peptide_df = pd.read_csv(peptide_file, sep='\t')
    protein_df = pd.read_csv(protein_file, sep='\t')

    result_columns = [
        'Protein', 'Num_peptides_both', 'Num_peptides_MSGF', 'Num_peptides_PA',
        'Num_unique_both', 'Num_unique_MSGF', 'Num_unique_PA',
        'Num_specs_both', 'Num_specs_MSGF', 'Num_specs_PA'
    ]
    result_df = pd.DataFrame(columns=result_columns)

    result_list = []
    for protein in protein_df['nextprot_accession']:
        matched_peptides = peptide_df[peptide_df['List_proteins'].str.contains(protein, na=False)]
        num_peptides_both = len(matched_peptides[matched_peptides['Num_specs_both'] > 0])
        num_peptides_MSGF = len(matched_peptides[matched_peptides['Num_specs_MSGF'] > 0])
        num_peptides_PA = len(matched_peptides[matched_peptides['Num_specs_PA'] > 0])
        num_unique_both = len(matched_peptides[(matched_peptides['Num_specs_both'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
        num_unique_MSGF = len(matched_peptides[(matched_peptides['Num_specs_MSGF'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
        num_unique_PA = len(matched_peptides[(matched_peptides['Num_specs_PA'] > 0) & (matched_peptides['Num_genes_saap'] == 1)])
        num_specs_both = matched_peptides['Num_specs_both'].sum()
        num_specs_MSGF = matched_peptides['Num_specs_MSGF'].sum()
        num_specs_PA = matched_peptides['Num_specs_PA'].sum()
        result_list.append({
            'Protein': protein,
            'Num_peptides_both': num_peptides_both,
            'Num_peptides_MSGF': num_peptides_MSGF,
            'Num_peptides_PA': num_peptides_PA,
            'Num_unique_both': num_unique_both,
            'Num_unique_MSGF': num_unique_MSGF,
            'Num_unique_PA': num_unique_PA,
            'Num_specs_both': num_specs_both,
            'Num_specs_MSGF': num_specs_MSGF,
            'Num_specs_PA': num_specs_PA
        })

    result_df = pd.concat([pd.DataFrame([result]) for result in result_list], ignore_index=True)
    result_df.to_csv(output_file, sep='\t', index=False)


def main():
    spectrum_file = 'example_reanalysis_spectrum.tsv'
    peptideatlas_file = 'PeptideAtlas_peptides.tsv'
    fasta_file = 'uniprotkb_human_proteome_UP000005640_with_isoforms_2024-10-08.fasta'
    output_file_peptide = 'peptide_out.tsv'
    print('files loaded')

    def get_peptide_sequence(row):
        return row['PeptideAtlas_peptide_demod'] if pd.notna(row['PeptideAtlas_peptide_demod']) else row['opt_global_UnmodPep']

    if not os.path.exists(output_file_peptide):
        with open(output_file_peptide, 'w') as f:
            f.write('')  # Create an empty file

    chunk_size = 1000
    total_chunks = sum(1 for _ in pd.read_csv(spectrum_file, sep='\t', chunksize=chunk_size))
    for i, chunk in enumerate(pd.read_csv(spectrum_file, sep='\t', chunksize=chunk_size)):
        if os.path.exists(output_file_peptide):
            if os.path.getsize(output_file_peptide) > 0:
                existing_output_df = pd.read_csv(output_file_peptide, sep='\t')
            else:
                existing_output_df = pd.DataFrame(columns=['Peptide sequence'])
            existing_peptides = set(existing_output_df['Peptide sequence'].unique())
        else:
            existing_peptides = set()

        chunk['Peptide sequence'] = chunk.apply(get_peptide_sequence, axis=1)
        chunk_peptides = chunk['Peptide sequence'].unique()
        unique_peptides = [pep for pep in chunk_peptides if pep not in existing_peptides]

        if unique_peptides:
            compare_peptide(spectrum_file, peptideatlas_file, fasta_file, output_file_peptide)
        
        print(f'Processed chunk {i + 1} of {total_chunks}')

if __name__ == '__main__':
    main()