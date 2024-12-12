# import json
# import re

# # Read the content of dataset_files.jsp
# with open('dataset_files.jsp', 'r') as file:
#     dataset_files = file.read()

# # Extract the JSON part from the dataset_files.jsp content
# json_data = re.search(r'var dataset_files = ({.*?});', dataset_files, re.DOTALL).group(1)

# # Load the JSON data
# data = json.loads(json_data)

# # Extract file names that end with .mzXML
# file_names = []
# for row in data['row_data']:
#     file_descriptor = row['file_descriptor']
#     match = re.search(r'ccms_peak/.*\.mzXML$', file_descriptor)
#     if match:
#         file_names.append(match.group())
#         print(file_names[-1])

# # Write the file names to mzXML_files.txt
# with open('mzXML_files.txt', 'w') as f:
#     for name in file_names:
#         f.write(name + '\n')

# Read the content of mzXML_files.txt

##################################################
# with open('mzml_files.txt', 'r') as file:
#     lines = file.readlines()

# # Remove 'ccms_peak/' from each line
# cleaned_lines = [line.replace('ccms_peak/', '') for line in lines]

# # Write the cleaned lines to mzml_files.txt
# with open('mzml_files.txt', 'w') as file:
#     file.writelines(cleaned_lines)

##################################################
with open('mzml_files.txt', 'r') as file:
    lines = file.readlines()

# Replace all '/' with '_'
cleaned_lines = [line.replace('/', '_') for line in lines]

# Write the cleaned lines to mzml_files.txt
with open('mzml_files.txt', 'w') as file:
    file.writelines(cleaned_lines)
