# Branchwater Test Data

This directory contains test datasets for PREPARE Branchwater functionality.

## Files

### `example.csv`
Original small synthetic test file with 2 records for basic testing.
- **Records:** 2
- **Organism:** Francisella adeliensis
- **Use case:** Unit tests, basic functionality validation

### `GCF_000006945.2_branchwater.csv` 
Real Branchwater results for **Salmonella enterica subsp. enterica** (RefSeq: GCF_000006945.2).
- **Records:** 8,027
- **Source:** Actual Branchwater analysis results
- **Coverage:** Worldwide samples from 114 countries
- **Containment range:** 0.100 - 1.000
- **Assay types:** 13 different types (WGS, OTHER, WGA, etc.)
- **Organisms:** 194 different organism types
- **Use case:** Performance testing, comprehensive analysis validation

### `salmonella_subset.csv`
Representative stratified sample from the full Salmonella dataset.
- **Records:** 83
- **Sampling strategy:** 
  - High containment (>0.8): 15 samples
  - Medium containment (0.3-0.8): 20 samples
  - Low containment (<0.3): 15 samples
  - Geographic diversity: Top 3 countries (China, USA, Malawi)
  - Assay type diversity: WGS, OTHER, WGA
  - Organism diversity: 26 different types
- **Use case:** Integration testing, moderate-scale analysis validation

### `salmonella_mini.csv`
Small representative sample for fast testing.
- **Records:** 19
- **Containment range:** 0.100 - 0.980
- **Countries:** 8 different countries
- **Assay types:** 3 types
- **Use case:** Quick tests, CI/CD pipelines, development

## Data Structure

All Branchwater CSV files follow this structure:

```csv
acc,containment,cANI,biosample,bioproject,assay_type,collection_date_sam,geo_loc_name_country_calc,organism,lat_lon
```

### Column Descriptions

- **acc:** SRA accession identifier
- **containment:** Containment value (0.0-1.0)
- **cANI:** Average nucleotide identity
- **biosample:** BioSample accession
- **bioproject:** BioProject accession  
- **assay_type:** Sequencing assay type
- **collection_date_sam:** Sample collection date
- **geo_loc_name_country_calc:** Calculated country name
- **organism:** Organism/sample type
- **lat_lon:** Geographic coordinates (when available)

## Usage Recommendations

- **Development/Unit Tests:** Use `salmonella_mini.csv` (19 records)
- **Integration Tests:** Use `salmonella_subset.csv` (83 records)  
- **Performance Tests:** Use `GCF_000006945.2_branchwater.csv` (8,027 records)
- **Basic Tests:** Use `example.csv` (2 records)

## Data Characteristics

### Geographic Distribution (Subset)
- Global representation across 12+ countries
- Major regions: USA, Europe, Asia, Africa
- Diverse sample sources and environments

### Containment Distribution (Subset)
- **High (>0.8):** ~18% of samples
- **Medium (0.3-0.8):** ~60% of samples  
- **Low (<0.3):** ~22% of samples

### Sample Types (Subset)
- Environmental metagenomes (wastewater, gut, food)
- Clinical isolates
- Synthetic/reference samples
- Various sequencing technologies

## Generation

The subset files were created using stratified sampling to ensure:
1. **Representative containment distribution**
2. **Geographic diversity**
3. **Assay type coverage**
4. **Organism type variety**
5. **Manageable file sizes for testing**

Random seed: 42 (subset), 123 (mini) for reproducible sampling.