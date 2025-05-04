import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import os

# Load the data
chunk_size = 10000
fasta_file = 'uniprotkb_human_proteome_UP000005640_with_isoforms_2024-10-08.fasta'

# Define a function to get the peptide sequence
def get_peptide_sequence(row):
    return row['PeptideAtlas_peptide_demod'] if pd.notna(row['PeptideAtlas_peptide_demod']) else row['sequence']

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

# Function to process a peptide
def process_peptide(to_parallel_process):
    peptide, peptide_data, peptideatlas_df = to_parallel_process
    
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
    #print(peptide_row)
    return peptide_row
    # Function to process an existing peptide
# def process_existing_peptide(peptide, peptide_data, peptideatlas_df):
#     # Build the output row without counting matches
#     peptide_row = {
#         'Peptide sequence': peptide,
#         'Peptide charge': peptide_data.apply(get_peptide_charge, axis=1).iloc[0],
#         'Protein identifier': peptide_data.apply(get_protein_id_from_msgf, axis=1).iloc[0],
#         'Num_specs_both': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
#         'Num_specs_MSGF': len(peptide_data[(pd.isna(peptide_data['PeptideAtlas_USI'])) & (pd.notna(peptide_data['sequence']))]),
#         'Num_specs_PA': len(peptide_data[(pd.notna(peptide_data['PeptideAtlas_USI'])) & (pd.isna(peptide_data['sequence']))]),
#         'PA_peptide': 1 if not peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide].empty else 0,
#         'PA_psms': len(peptideatlas_df[peptideatlas_df.iloc[:, 5] == peptide]),
#         'Num_proteins': None,
#         'List_proteins': None,
#         'Num_genes': None,
#         'List_genes': None,
#         'Num_proteins_saap': None,
#         'List_proteins_saap': None,
#         'Num_genes_saap': None,
#         'List_genes_saap': None
#     }
#     print(f"Existing Peptide ID: {peptide_row['Peptide sequence']}")
#     return peptide_row
# Load the PeptideAtlas data
peptideatlas_df = pd.read_csv('PeptideAtlas_peptides.tsv', sep='\t')

# Use a set to keep track of unique peptides
unique_peptides_set = set()

if __name__ == "__main__":
    # Open the output file
    peptide_file_exists = os.path.exists('example_reanalysis_peptide.tsv')
    with open('example_reanalysis_peptide.tsv', 'a+') as output_file:
        # Write the header
        if not peptide_file_exists:
            output_file.write('\t'.join([
                'Peptide sequence', 'Peptide charge', 'Protein identifier', 'Num_specs_both', 'Num_specs_MSGF', 'Num_specs_PA',
                'PA_peptide', 'PA_psms', 'Num_proteins', 'List_proteins', 'Num_genes', 'List_genes', 'Num_proteins_saap',
                'List_proteins_saap', 'Num_genes_saap', 'List_genes_saap'
            ]) + '\n')
        
        # Read the spectrum file in chunks
        # Check if the file exists
        if peptide_file_exists:
            # Read the existing peptides from the file
            with open('example_reanalysis_peptide.tsv', 'r') as existing_file:
                # Skip the header
                next(existing_file)
                for line in existing_file:
                    peptide = line.split('\t')[0]  # Extract the peptide sequence (first column)
                    unique_peptides_set.add(peptide)
        print(f"Number of unique peptides added: {len(unique_peptides_set)}")
        for chunk in pd.read_csv('example_reanalysis_spectrum.tsv', sep='\t', chunksize=chunk_size):
            chunk['Peptide sequence'] = chunk.apply(get_peptide_sequence, axis=1)
            counter = 0
            total_peptides = len(chunk['Peptide sequence'].unique())

            def update_counter(result):
                global counter
                counter += 1
                # print(f"Processed {counter}/{total_peptides} peptides")

            with Pool(cpu_count()) as pool:
                peptides_to_process = [
                    (peptide, chunk[chunk['Peptide sequence'] == peptide], peptideatlas_df)
                    for peptide in chunk['Peptide sequence'].unique()
                ]
                
                results = []
                to_parallel_process = []
                # not_parallel_process = []

                for peptide, peptide_data, pa_df in peptides_to_process:
                    if peptide not in unique_peptides_set:
                        #print(f"Appending New Peptide: {peptide}")
                        unique_peptides_set.add(peptide)
                        to_parallel_process.append((peptide, peptide_data, pa_df))
                    else:
                        #print(f"Processing Existing Peptide: {peptide}")
                        # not_parallel_process.append((peptide, peptide_data, pa_df))
                        pass

                # Parallel process the new peptides
                if to_parallel_process:
                    for result in pool.imap_unordered(process_peptide, to_parallel_process):
                        results.append(result)
                        update_counter(result)
                        print(f"Processing new peptide: {result['Peptide sequence']} ({counter}/{len(to_parallel_process)}) in chunk {chunk.index[0] // chunk_size + 1}")
                        # Calculate and print progress percentage
                        total_new_peptides = len(to_parallel_process)
                        total_peptides_all_chunks = len(unique_peptides_set)
                        percentage = (counter / total_new_peptides) * 100 if total_new_peptides > 0 else 100
                        overall_percentage = (len(unique_peptides_set) / (total_peptides_all_chunks + total_new_peptides - counter)) * 100
                        print(f"Progress: {percentage:.2f}% of new peptides in current chunk, {overall_percentage:.2f}% of all peptides processed so far")

                # Process the existing peptides
                print('Processing existing peptides...SIKE')
                # for peptide, peptide_data, pa_df in not_parallel_process:
                #     results.append(process_existing_peptide(peptide, peptide_data, pa_df))
                #     update_counter(None)
                    #print(f"Processed {counter - len(not_parallel_process)}/{len(not_parallel_process)} old peptides in chunk {chunk.index[0] // chunk_size + 1}")
                    
                for peptide_row in results:
                    # Write the results to the output file
                    # Write the peptide row as a tab-separated line
                    row_values = [
                        str(peptide_row.get(col, '')) for col in [
                            'Peptide sequence', 'Peptide charge', 'Protein identifier', 'Num_specs_both', 'Num_specs_MSGF', 'Num_specs_PA',
                            'PA_peptide', 'PA_psms', 'Num_proteins', 'List_proteins', 'Num_genes', 'List_genes', 'Num_proteins_saap',
                            'List_proteins_saap', 'Num_genes_saap', 'List_genes_saap'
                        ]
                    ]
                    output_file.write('\t'.join(row_values) + '\n')
                    output_file.flush()  # Ensure the data is written immediately
                    

                pool.close()
                pool.join()

            print(f'Processed chunk {chunk.index[0] // chunk_size + 1}')
#128 2.06am  8m 50 peptides ->  16m 100 peptieds -> 16*120 min 12000 peptides
#178 2.14am
#228 2.19am 
#328 2.31am