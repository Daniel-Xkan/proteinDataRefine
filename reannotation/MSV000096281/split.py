import csv

# Read the input file
with open('mzml_files.txt', 'r') as infile:
    lines = infile.readlines()

# Split each line by '_'
split_lines = [line.strip().split('_') for line in lines]

# Write the split data to a CSV file
with open('split.csv', 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerows(split_lines)