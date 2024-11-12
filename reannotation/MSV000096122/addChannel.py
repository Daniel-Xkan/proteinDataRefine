import pandas as pd

# Read the CSV files
set_label_df = pd.read_csv('setLable.csv')
table_source_df = pd.read_csv('tableSource.csv')

# Create a new DataFrame to store the modified rows
new_rows = []

# Iterate over each row in tableSource.csv
for _, row in table_source_df.iterrows():
    experiment = row['Experiment']
    set_number = int(experiment.split('Set')[-1])
    
    # Filter setLable.csv to find matching TMT Batch Number
    matching_set_label = set_label_df[set_label_df['TMT Batch Number'].str.contains(f'TMT{set_number}')]
    
    if not matching_set_label.empty:
        # Get the number of duplicates needed
        num_duplicates = len(matching_set_label)
        
        for i in range(num_duplicates):
            new_row = row.copy()
            new_row['Channel'] = matching_set_label.iloc[i]['TMT Channel']
            new_row['BioReplicate'] = matching_set_label.iloc[i]['Donors']
            new_rows.append(new_row)
        
        # Add the 'Reference' row
        reference_row = row.copy()
        reference_row['Channel'] = 131
        reference_row['BioReplicate'] = 'Reference'
        new_rows.append(reference_row)

# Create a DataFrame from the new rows
new_table_source_df = pd.DataFrame(new_rows)

# Output the altered DataFrame to a new CSV file
new_table_source_df.to_csv('MSV000096122_annotated.csv', index=False)