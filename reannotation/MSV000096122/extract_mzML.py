import json
import re

# Read the JSP file
with open('dataset_files.jsp', 'r') as file:
    jsp_content = file.read()

# Attempt to match the dataset_files variable with a more flexible regex pattern
json_match = re.search(r'dataset_files\s*=\s*({.*?})\s*;\s*', jsp_content, re.DOTALL)

# Check if the match was successful
if json_match:
    # Extract and clean up JSON-like content
    json_content = json_match.group(1)
    json_content = re.sub(r',\s*}', '}', json_content)  # Remove trailing commas before closing braces
    json_content = re.sub(r',\s*\]', ']', json_content)  # Remove trailing commas before closing brackets
    
    try:
        # Parse the JSON content
        data = json.loads(json_content)
        
        # Extract all .mzML file names
        mzml_files = [item['name'] for item in data['row_data'] if item['name'].endswith('.mzML')]
        
        # Write the filenames to a text file
        with open('mzml_files.txt', 'w') as output_file:
            for filename in mzml_files:
                output_file.write(filename + '\n')
        
        print("mzML filenames successfully written to mzml_files.txt")
    except json.JSONDecodeError as e:
        print("JSON decoding error:", e)
else:
    # Print the JSP content to help with debugging
    print("No JSON data found in JSP content. Check the structure of the file.")