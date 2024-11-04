import pandas as pd
import re

# Load filenames
with open('mzml_files.txt', 'r') as file:
    filenames = file.read().splitlines()

# Define the structure of the DataFrame with the required columns
columns = ["Dataset", "Run", "Cohort", "Subject", "Condition", "TimePoint", "BioReplicate", 
           "Experiment", "Batch", "Channel", "Fraction", "Species", "Disease", "Tissue", "Enzyme", "Notes"]
data = {col: [] for col in columns}

# Pattern for filenames with more identifiable structure
pattern = re.compile(
    r'(?P<Dataset>\d{8})_'           # Date (Dataset)
    r'(?P<Run>QEh\d+_\w+)_?'         # Run (e.g., QEh1_LC1)
    r'(?P<Cohort>\w+)?_'             # Cohort
    r'(?P<Subject>\w+)?_'            # Subject
    r'(?P<Condition>\w+)?_'          # Condition
    r'(?P<TimePoint>\d+)?_'          # TimePoint
    r'(?P<BioReplicate>R\d+)?'       # BioReplicate (e.g., R1, R2)
    r'\.mzML'                        # File extension
)

# Parse each filename and populate columns
for filename in filenames:
    match = pattern.match(filename)
    if match:
        # Fill matched parts
        for col, value in match.groupdict().items():
            data[col].append(value if value else "")
        # Fill remaining columns with empty strings
        for col in columns:
            if col not in match.groupdict():
                data[col].append("")
    else:
        # If the filename doesn't match the pattern, add blanks
        for col in columns:
            data[col].append("")

# Create DataFrame
df = pd.DataFrame(data)

# Save to TSV
df.to_csv("reannotated_MSV000096130.tsv", sep="\t", index=False)

print("Parsed TSV saved as parsed_mzml_files.tsv")
