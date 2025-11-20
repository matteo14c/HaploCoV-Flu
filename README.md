# HaploCoV -- HC_comb Grouping and Amino-Acid Annotation for IAV 2.3.4.4b

Scripts and Pipeline for Chiara et al. (submitted)

This repository contains all scripts used to generate **HC_comb groups**
and the **annotation of amino-acid substitutions** across the **eight
genomic segments** of *Influenza A Virus (IAV) clade 2.3.4.4b*, as
performed in **Chiara et al. (submitted)**.

Because of GISAID restrictions, **the raw sequences used in the analysis
cannot be distributed directly**.\
Users must obtain these sequences independently using the **EPISET
referenced in the manuscript**.

Once the input sequences are obtained, the provided pipeline
(`GenerateResults.sh`) reproduces all tables and annotations described
in the manuscript.

------------------------------------------------------------------------

## 1. Input Data Requirements

GISAID sequences must be retrieved using the EPISET referenced in the
paper.\
All sequences must be placed into a **single multi-FASTA file** named:

    fasta/Sequences.fasta

### **FASTA header format**

Each sequence header **must** follow this structure:

    >A_fox_England_015850_2022|EPI_ISL_17072388|NP

Where:

-   **virus name** (e.g., `A_fox_England_015850_2022`)
-   followed by `|`
-   **GISAID EPI ID** (e.g., `EPI_ISL_17072388`)
-   followed by `|`
-   **segment name** (`PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `M`, or
    `NS`)

This format is required for correct parsing and grouping by the
pipeline.

------------------------------------------------------------------------

## 2. Running the Pipeline

Once `Sequences.fasta` is placed in `fasta/`, run:

``` bash
bash GenerateResults.sh
```

This script regenerates **all results** described in the manuscript.

------------------------------------------------------------------------

## 3. Output Files

The pipeline produces:

### **1. HC_comb_file.tsv**

-   HaploCoV grouping results for each sequence\
-   Corresponds to **Supplementary Table 1** in the manuscript\
-   Lists the *HC_comb designations* for all isolates

### **2. Eight annotated phenetic-pattern files**

One per segment, named:

    PB2.annot_V2.csv
    PB1.annot_V2.csv
    PA.annot_V2.csv
    HA.annot_V2.csv
    NP.annot_V2.csv
    NA.annot_V2.csv
    M.annot_V2.csv
    NS.annot_V2.csv

Each file contains:\
- phenetic patterns (presence/absence of nonsynonymous substitutions)\
- annotated amino-acid effects\
- These correspond to **Supplementary Table S6** in the manuscript

------------------------------------------------------------------------

## 4. Software Requirements

### **Required external tool**

HaploCoV relies on:

-   **nucmer** (from the MUMmer suite)

Please install it before running the pipeline.\
For example on Debian/Ubuntu:

``` bash
sudo apt-get install mummer
```

Or from source:

https://mummer4.github.io/

### **Required language**

A working **Perl interpreter** is needed, since all HaploCoV scripts are
written in Perl.

------------------------------------------------------------------------

## 5. Repository Structure

    scripts/            # All Perl scripts used in the analysis
    fasta/              # Contains Sequences.fasta provided by user
    data/               # Intermediate and output files
    reference/          # Reference segment FASTA files
    annotationFiles/    # Segment-specific annotation tables
    GenerateResults.sh  # Master pipeline script
    README.md           # This file

------------------------------------------------------------------------

## 6. Citation

If you use this code, please cite:

**Chiara et al., (submitted).**\
"Patterns of evolutionary change in IAV clade 2.3.4.4b revealed by
HaploCoV."

------------------------------------------------------------------------

For questions, comments, or issues, please open a GitHub issue.\
We welcome community feedback.
