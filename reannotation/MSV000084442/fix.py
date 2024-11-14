import pandas as pd

# # Load the annotated CSV file
# df = pd.read_csv('MSV000084442_annotated.csv')

# # Define the mapping from Clinical ID to Sample ID and Cancer type
# mapping = {
#     'DFCI-5341': ('CLL A', 'CLL'),
#     'DFCI-5328': ('CLL B', 'CLL'),
#     'DFCI-5283': ('CLL C', 'CLL'),
#     '13240-002': ('MEL 1', 'Melanoma'),
#     '13240-005': ('MEL2', 'Melanoma'),
#     '13240-006': ('MEL3', 'Melanoma'),
#     '13240-015': ('MEL15', 'Melanoma'),
#     'CP-594_v1': ('OV1', 'Ovarian'),
#     '14362-007': ('GBM7', 'GBM'),
#     'H4198 BT187': ('GBM9', 'GBM'),
#     'H4512 BT145': ('GBM11', 'GBM'),
#     'Pat9': ('', 'ccRCC')
# }

# # Update the Bioreplicate and Cohort columns based on the Run column
# for clinical_id, (sample_id, cancer_type) in mapping.items():
#     df.loc[df['Run'].str.contains(clinical_id, na=False), 'Bioreplicate'] = sample_id
#     df.loc[df['Run'].str.contains(clinical_id, na=False), 'Cohort'] = cancer_type

# # Save the updated CSV file
# df.to_csv('MSV000084442_annotated.csv', index=False)

# Load the annotated CSV file

df = pd.read_csv('MSV000084442_annotated.csv')

# Update the Experiment column based on the Run column
df['Experiment'] = df['Run'].apply(lambda x: x.split('/')[2] if '/' in x else '')

# Save the updated CSV file
df.to_csv('MSV000084442_annotated.csv', index=False)