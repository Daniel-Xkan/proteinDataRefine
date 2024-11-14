import csv
import os
import re

# Read mzML file names
with open('mzml_files.txt', 'r') as f:
    mzml_files = f.read().splitlines()

# Create separated.csv
with open('separated.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for file in mzml_files:
        parts = file.split('/')
        filename = parts[-1]
        elements = filename.split('_')
        writer.writerow(parts[:-1] + elements)

# Read cohort data
cohort_data = {
    "CLL A": ["A0301", "A3101", "B1402", "B3502", "C0401", "C0802"],
    "CLL B": ["A0206", "A2402", "B0801", "B5101", "C0702", "C1402"],
    "CLL C": ["A0101", "A0201", "B0702", "B0801", "C0701", "C0702"],
    "GBM7": ["A2402", "A6801", "B3502", "B3503", "C0401", "C0401"],
    "GBM7 +IFN": ["A2402", "A6801", "B3502", "B3503", "C0401", "C0401"],
    "GBM9": ["A0101", "A3201", "B1501", "B1501", "C0303", "C0304"],
    "GBM9 +IFN": ["A0101", "A3201", "B1501", "B1501", "C0303", "C0304"],
    "GBM11": ["A0101", "A3201", "B1302", "B4002", "C0202", "C0602"],
    "GBM11 +IFN": ["A0101", "A3201", "B1302", "B4002", "C0202", "C0602"],
    "MEL 1": ["A0201", "A2402", "B1501", "B4402", "C0501", "C0702"],
    "MEL 2": ["A0101", "A0101", "B3801", "B5601", "C0102", "C0602"],
    "MEL 2 +IFN": ["A0101", "A0101", "B3801", "B5601", "C0102", "C0602"],
    "MEL 3": ["A0201", "A0301", "B2705", "B4701", "C0102", "C0602"],
    "MEL 3 +IFN": ["A0201", "A0301", "B2705", "B4701", "C0102", "C0602"],
    "MEL 15": ["A0201", "A0202", "B1302", "B4002", "C0202", "C0602"],
    "MEL 15 +IFN": ["A0201", "A0202", "B1302", "B4002", "C0202", "C0602"],
    "OV1": ["A0201", "A2402", "B3503", "B4402", "C0501", "C1203"],
    "Mel8": ["A0101", "A0301", "B0702", "B0801", "C0701", "C0702"],
    "Mel12": ["A0101", "B0801", "C0701"],
    "Mel15": ["A0301", "A6801", "B2705", "B3503", "C0202", "C0401"],
    "Mel16": ["A0101", "A2402", "B0702", "B0801", "C0701", "C0702"],
    "OvCa9": ["A0201", "A0301", "B0702", "B4002", "C0702", "C1201"],
    "OvCa10": ["A0201", "A1101", "B4405", "B5101", "C0202", "C1502"],
    "OvCa12": ["A2402", "A3101", "B3503", "B4901", "C0701", "C1203"],
    "OvCa15": ["A1101", "A2402", "B0702", "B5501", "C0303", "C0702"],
    "OvCa28": ["A0101", "A0201", "B2705", "B5201", "C0102", "C0202"],
    "OvCa39": ["A2501", "A3101", "B0702", "B1801", "C1203", "C0702"],
    "OvCa48": ["A0201", "A2501", "B1501", "B4102", "C0304", "C1701"],
    "OvCa54": ["A0201", "A1101", "B3501", "B3503", "C0401", "C1203"],
    "OvCa60": ["A2402", "A2501", "B1302", "B1801", "C1203", "C0602"],
    "OvCa66": ["A1101", "A2902", "B1801", "B4403", "C0501", "C1601"],
    "OvCa68": ["A0201", "A0101", "B4402", "B3701", "C0602", "C0501"],
    "OvCa72": ["A0301", "A0101", "B0801", "B0702", "C0702", "C0701"],
    "OvCa74": ["A0201", "B1801", "B5101", "C0702", "C1502"],
    "OvCa79": ["A0101", "A3101", "B0801", "B5101", "C0701", "C1502"],
    "OvCa80": ["A2501", "A3201", "B1801", "B3901", "C1203"],
    "OvCa81": ["A0201", "B4501", "B5601", "C0702", "C0102"],
    "OvCa82": ["A0101", "A0301", "B0801", "B3801", "C0701", "C1203"],
    "OvCa84": ["A0201", "B0702", "B4402", "C0702", "C0501"],
    "OvCa99": ["A0201", "A2402", "B1302", "B4001", "C0304", "C0602"],
    "OvCa100": ["A0201", "B0702", "B4102", "C0702", "C1701"],
    "OvCa103": ["A0201", "A2402", "B2702", "B2705", "C0202"],
    "OvCa104": ["A0301", "B0702", "B3508", "C0401", "C0702"],
    "OvCa105": ["A2601", "A6801", "B1801", "B5501", "C0303", "C0701"],
    "OvCa109": ["A0201", "A2301", "B4001", "B4901", "C0701", "C0304"],
    "OvCa111": ["A0101", "A2501", "B0801", "B4402", "C0501", "C0701"],
    "OvCa114": ["A2902", "B4403", "C1601"],
    "OvCa73": ["A0101", "B0801", "C0701"]
}

# Create MSV000084172_annotated.csv
with open('separated.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    with open('MSV000084172_annotated.csv', 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Dataset', 'Run', 'Cohort', 'Subject', 'Condition', 'TimePoint', 'BioReplicate', 'Experiment', 'Batch', 'Channel', 'Fraction', 'Species', 'Disease', 'Tissue', 'Enzyme', 'Notes'])
        
        for row in reader:
            dataset = "MSV000084172"
            run = '/'.join(row[:3]) + '/' + '_'.join(row[3:])
            timepoint = row[1]
            experiment = row[-3]
            batch = row[-2] if len(row) > 4 else ''
            
            # Determine Cohort and BioReplicate
            cohort = ''
            bioreplicate = ''
            pattern = re.compile(r'[A-C]\d{4}')
            for key, values in cohort_data.items():
                if any(value in run for value in values):
                    cohort = key
                    bioreplicate = values[0]  # Assuming first value is BioReplicate
                    break
            else:
                match = pattern.search(run)
                if match:
                    cohort = match.group(0)
            
            # Write to annotated CSV
            writer.writerow([dataset, run, cohort, cohort, '', timepoint, bioreplicate, experiment, batch, '', '', '', '', '', '', ''])

# Open the file again to update Disease, Batch, and Experiment
with open('MSV000084172_annotated.csv', 'r') as infile:
    reader = csv.reader(infile)
    rows = list(reader)

header = rows[0]
data_rows = rows[1:]

# Define the mapping for Disease, Sample ID, and Experiment
cancer_mapping = {
    "CLL A": ("CLL", "DFCI-5341", "CLL, CD19 sorted"),
    "CLL B": ("CLL", "DFCI-5328", "CLL, CD19 sorted"),
    "CLL C": ("CLL", "DFCI-5283", "CD19 sorted"),
    "MEL 1": ("Melanoma", "13240-002", "cells"),
    "MEL 2": ("Melanoma", "13240-005", "cells"),
    "MEL 2 +IFN": ("Melanoma", "13240-005", "cells + IFNg"),
    "MEL 3": ("Melanoma", "13240-006", "cells"),
    "MEL 3 +IFN": ("Melanoma", "13240-006", "cells + IFNg"),
    "MEL 15": ("Melanoma", "13240-015", "cells"),
    "MEL 15 +IFN": ("Melanoma", "13240-015", "cells + IFNg"),
    "OV1": ("Ovarian", "CP-594_v1", "cells"),
    "GBM7": ("GBM", "14362-007", "cells"),
    "GBM7 +IFN": ("GBM", "14362-007", "cells + IFNg"),
    "GBM9": ("GBM", "H4198 BT187", "cells"),
    "GBM9 +IFN": ("GBM", "H4198 BT187", "cells+IFNg"),
    "GBM11": ("GBM", "H4512 BT145", "cells"),
    "GBM11 +IFN": ("GBM", "H4512 BT145", "cells + IFNg"),
    "Pat9": ("ccRCC", "Pat9", "cells"),
    "Pat9 +IFN": ("ccRCC", "Pat9", "cells+ IFNg")
}

# Update the rows with the mapping
for row in data_rows:
    cohort = row[2]
    if cohort in cancer_mapping:
        disease, sample_id, experiment = cancer_mapping[cohort]
        row[12] = disease
        row[8] = sample_id
        row[7] = experiment

# Write the updated rows back to the file
with open('MSV000084172_annotated.csv', 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(header)
    writer.writerows(data_rows)
