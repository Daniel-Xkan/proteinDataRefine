import pandas as pd
import re

# Read file names from mzml_files.txt
with open("mzml_files.txt", "r") as file:
    file_names = [line.strip() for line in file.readlines()]

# Initialize lists to hold extracted components
datasets = []
cohorts = []
experiments = []
batches = []
fractions = []

# Regular expression pattern to parse the filename
# It covers cases with and without "rern" in the set name
pattern = re.compile(r"(ScltlMscl)(Set\d+)(rern)?(BRPhs)(Fr\d+)\.mzML")

# Parse each filename and extract components
for filename in file_names:
    match = pattern.match(filename)
    if match:
        datasets.append("MSV000096122")       # Dataset as provided in interpretation
        cohorts.append(match.group(1))        # Cohort: "ScltlMscl" (skeletal mass spec-based analysis)
        experiments.append(match.group(2))    # Experiment: Set number
        batch_prefix = "rern" if match.group(3) else ""
        batches.append(batch_prefix + match.group(4))  # Batch: either "BRPhs" or "rernBRPhs"
        fractions.append(match.group(5))      # Fraction
    else:
        # If format does not match, add empty strings for consistency
        datasets.append("")
        cohorts.append("")
        experiments.append("")
        batches.append("")
        fractions.append("")

# Create DataFrame with structured data
df = pd.DataFrame({
    "Dataset": datasets,
    "Run": file_names,
    "Cohort": cohorts,
    "Subject": ["" for _ in file_names],  # Placeholder, as subject info isn't available
    "Condition": ["" for _ in file_names],  # Placeholder
    "TimePoint": ["" for _ in file_names],  # Placeholder
    "BioReplicate": cohorts,  # Assuming BioReplicate aligns with cohort info based on the interpretations
    "Experiment": experiments,
    "Batch": batches,
    "Channel": ["" for _ in file_names],  # Placeholder
    "Fraction": fractions,
    "Species": ["" for _ in file_names],  # Placeholder
    "Disease": ["" for _ in file_names],  # Placeholder
    "Tissue": ["" for _ in file_names],  # Placeholder
    "Enzyme": ["" for _ in file_names],  # Placeholder
    "Notes": ["" for _ in file_names],  # Placeholder
})

# Save the structured data to a CSV file with empty values for missing fields
df.to_csv("MSV000096122_ skeletal_muscle_annotation.csv", sep=',', index=False)

