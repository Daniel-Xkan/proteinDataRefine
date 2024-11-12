import re
import csv
def parse_filename(filename):
    pattern = re.compile(r'(\d{6})_(\w+)_([A-Z0-9-]+)_([A-Za-z]+)_([A-Za-z0-9-]+)_([0-9]+)_([A-Za-z0-9_]+)_([0-9]+[a-z]*)_([0-9-]+mz)_msms([0-9A-Za-z]+)_standard\.mzML(?:\.scans)?')
    match = pattern.match(filename)
    if match:
        date, _, bio_replicate, tissue, batch, _, experiment, experiment2, experiment3, fraction = match.groups()
        return {
            'Dataset': 'MSV000096271',
            'Run': filename,
            'Cohort': tissue,
            'Subject': tissue,
            'Condition': '',
            'TimePoint': date,
            'BioReplicate': bio_replicate,
            'Experiment': experiment+'_'+experiment2+'_'+experiment3,
            'Batch': batch,
            'Channel': '',
            'Fraction': 'msms'+fraction,
            'Species': '',
            'Disease': '',
            'Tissue': tissue,
            'Enzyme': '',
            'Notes': ''
        }
    return None


def read_filenames(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

def write_csv(data, output_file):
    fieldnames = ['Dataset', 'Run', 'Cohort', 'Subject', 'Condition', 'TimePoint', 'BioReplicate', 'Experiment', 'Batch', 'Channel', 'Fraction', 'Species', 'Disease', 'Tissue', 'Enzyme', 'Notes']
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

def main():
    input_file = 'mzml_files.txt'
    output_file = 'MSV000096271_annotated.csv'
    filenames = read_filenames(input_file)
    data = []
    for filename in filenames:
        parsed_data = parse_filename(filename)
        if parsed_data:
            data.append(parsed_data)
    write_csv(data, output_file)

if __name__ == "__main__":
    main()