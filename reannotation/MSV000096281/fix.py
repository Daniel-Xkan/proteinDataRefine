import pandas as pd

# Read mzml_files.txt and split by '_', then save to split.csv
with open('mzml_files.txt', 'r') as file:
    lines = file.readlines()

split_data = [line.strip().split('_') for line in lines]
split_df = pd.DataFrame(split_data)
split_df.to_csv('split.csv', index=False, header=False)

# Load the original annotated CSV file
df = pd.read_csv('MSV000096281_annotated.csv')

# Function to generate the Fraction value
def generate_fraction(run):
    parts = run.split('_')
    if parts[-1].startswith(('01', '02', '03', '04', '05')):
        return f"{parts[-2]}_{parts[-1].replace('.mzML', '')}"
    return None

# Apply the function to each row and update the 'Fraction' column
df['Fraction'] = df['Run'].apply(generate_fraction)

# Save the updated DataFrame to the original CSV file
df.to_csv('MSV000096281_annotated.csv', index=False)