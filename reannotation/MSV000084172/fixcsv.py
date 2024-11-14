import csv
import re

# input_file = 'MSV000084172_annotated.csv'
# output_file = 'MSV000084172_annotated_fixed.csv'

# pattern = re.compile(r'[A-Z]\d{4}')

# with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
#     reader = csv.DictReader(infile)
#     fieldnames = reader.fieldnames
#     writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    
#     writer.writeheader()
#     for row in reader:
#         run_value = row['Run']
#         match = pattern.search(run_value)
#         if match:
#             row['BioReplicate'] = match.group(0)
#         writer.writerow(row)

input_file = 'MSV000084172_annotated.csv'
output_file = 'MSV000084172_annotated_fixed.csv'

pattern = re.compile(r'R\d_\d\d')

with open(input_file, mode='r', newline='') as infile, open(output_file, mode='w', newline='') as outfile:
    reader = csv.DictReader(infile)
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in reader:
        run_value = row['Run']
        match = pattern.search(run_value)
        if match:
            row['Fraction'] = match.group(0)
        writer.writerow(row)