# Automation for Protein and Spectrum Analysis

## 1. PeptideAtlas.org Coverage Extraction for Missing MassIVE-KB Proteins (/trial1) finished

### Objective
Identify peptides from PeptideAtlas.org not found in MassIVE-KB to improve protein coverage analysis.

### Steps
1. **Filter Proteins**:
   - Start with PeptideAtlas proteins (2024-01 build).
   - Skip proteins present in `MassIVE-KB_HPP_proteins.tsv`.
   - Access the PeptideAtlas page for remaining proteins via the canonical link.

2. **Extract Peptide Data**:
   - Parse "Distinct Observed Peptides" table.
   - Skip peptides present in `MassIVE-KB_HPP_peptides.tsv` or if `is_unique` is "N."
   - Save relevant peptide data, including:
     - Protein "Gene Accession."
     - Peptide details from "Distinct Observed Peptides."

3. **Parse Experiment Observations**:
   - Access the peptide page through the "Accession" link.
   - Parse "Observed in Experiments" table.
   - Output relevant data to `PeptideAtlas_peptide_observations.tsv`.
   - Track experiment count (`NumExperiments`) and save peptide data to `PeptideAtlas_peptides.tsv`.

## 2. Spectra Extraction for MassIVE Reanalysis (/trial2)

### Part 1: Import Datasets (inprogress)
1. **Dataset Preparation**:
   - Use PRIDE-IMPORT and dataset conversion tools.
   - Create annotation files per dataset, focusing on fields like `Dataset`, `Run`, `Cohort`, `BioReplicate`, and `Experiment`.

2. **Meta-workflow Execution**:
   - Run searches and extract relevant data to enhance the knowledge base (KB).

### Part 2: Peptide Spectrum Matches (/trial2/multithread) finished
1. **Prepare Spectrum Files**:
   - Create text files with PSM Universal Spectrum Identifiers (USIs) for peptides.
   - Generate MGF files containing spectra.

2. **Categorize Peptides**:
   - Generate `<protein_id>_unmodified_peptides_usi.txt` for unmodified peptides.
   - Generate `<protein_id>_modified_peptides_usi.txt` for modified peptides.

3. **Extract and Validate USIs**:
   - Visit PeptideAtlas peptide pages and extract USIs.
   - Discard invalid USIs and categorize based on peptide modifications.

4. **Exact-match Searches**:
   - Run searches against predicted spectra (Prosit) and MassIVE-KB to identify confirming or competing matches.
