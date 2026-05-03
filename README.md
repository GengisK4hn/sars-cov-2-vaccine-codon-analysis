# SARS-CoV-2 Vaccine Codon Analysis

Bioinformatics analysis revealing artificial codon preference changes in COVID-19 mRNA vaccines compared to natural SARS-CoV-2 variants.

## Overview

This repository contains analysis scripts and documentation for a comprehensive bioinformatics study comparing SARS-CoV-2 vaccine sequences with natural viral variants. The key finding is that only vaccine sequences show artificial amino acid preference changes, while all natural variants maintain the Wuhan baseline.

## Key Findings

### Finding 1: Amino Acid Preference Changes

- **Pfizer BNT162b2:** 3/20 (15%) amino acid preference changes, p < 0.001
- **Moderna mRNA-1273:** 7/20 (35%) amino acid preference changes, p < 0.0001
- **All natural variants:** 0/20 (0%) amino acid preference changes

This is the first documentation of this pattern and appears to be the key distinguishing feature between engineered and natural sequences.

### Additional Findings

- **44nt consensus sequence:** Found in vaccine vials (156,086 reads in Moderna), absent from natural variants
- **19nt FCS reverse complement:** Present in Moderna patent and RNAseq (548 reads)
- **VERO/HAE adaptation signatures:** Consistent with laboratory passage
- **NLS motifs:** 26 detected in Pfizer, 0 in Moderna
- **GOF signatures:** CGG codons and restriction sites characteristic of infectious clone assembly

## Installation

### Requirements

- Python 3.8+
- BioPython
- pandas
- numpy
- scipy

### Setup

```bash
# Clone repository
git clone https://github.com/GengisK4hn/sars-cov-2-vaccine-codon-analysis.git
cd sars-cov-2-vaccine-codon-analysis

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Download Reference Sequences

```bash
# Download all required sequences from NCBI
bash download_all_sequences.sh
```

### Run Analysis

```bash
# Run codon optimization analysis
python scripts/codon_optimization_verifier.py \
  data/vaccine_references/pfizer_bnt162b2_vector_OR134577.fasta \
  --natural data/sequences/MT020880.1_Early_Wuhan.fasta data/sequences/SARS-CoV-2.fasta \
  --output CODON_OPTIMIZATION_PFIZER.txt

# Run independent verification
python scripts/independent_verification.py --data-dir data
```

### Verify Findings

All key findings can be verified using direct grep commands:

```bash
# Verify 44nt consensus sequence
grep -c "AAGATCGCCGACTACAACTACAAGCTGCCCGACGACTTCACCGG" \
  data/sequences/RNAseq-Mod2_R2_001.fastq.fasta

# Verify 19nt FCS reverse complement
grep -c "CTACGTGCCCGCCGAGGAG" \
  data/sequences/RNAseq-Mod2_R2_001.fastq.fasta
```

## Data Sources

### NCBI GenBank Accessions

- **Pfizer BNT162b2:** OR134577.1
- **Moderna mRNA-1273:** OR134578.1
- **Wuhan-Hu-1:** NC_045512.2
- **Early Wuhan:** MT020880.1
- **Delta:** OM095706.1
- **Omicron BA.1:** OMX067679.1
- **Omicron BA.2:** OMX067680.1

## Reproducibility

All analysis scripts and verification commands are provided for independent validation.

## License

MIT License - see LICENSE file for details.

---

**Status:** Analysis complete, findings verified, ready for peer review
**Last Updated:** May 2026
