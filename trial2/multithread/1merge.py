import os

# Directory containing the files
directory = os.path.dirname(__file__)

# Output file
output_file = os.path.join(directory, 'zmerged.txt')

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file ends with 'usi.txt'
        if filename.endswith('usi.txt'):
            file_path = os.path.join(directory, filename)
            # Open each file in read mode
            with open(file_path, 'r') as infile:
                # Write the contents of the file to the output file
                outfile.write(infile.read())
                outfile.write('\n')  # Add a newline character after each file's content``