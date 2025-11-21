#!/usr/bin/env bash
set -euo pipefail

##############################
# Configuration
##############################

SEGMENTS=(PB2 PB1 PA HA NP NA M NS)

FASTA_DIR="fasta"
REF_DIR="reference"
DATA_DIR="data"
ANNOT_DIR="annotationFiles"
SCRIPTS_DIR="scripts"

MULTIFASTA="$FASTA_DIR/Sequences.fasta"
META_TABLE="$DATA_DIR/MetadataTable.tsv"

OUTPUT_TABLE="HC_comb_file.tsv"

# Default number of threads
THREADS=8

########################################
# Argument parsing
########################################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads|-t)
            THREADS="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--threads N]"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done


##############################
# Helper Functions
##############################

info()   { printf "\n[INFO] %s\n" "$1"; }
error()  { printf "\n[ERROR] %s\n" "$1" >&2; exit 1; }
banner() { printf "\n========== %s ==========\n" "$1"; }

check_file() {
    [[ -f "$1" ]] || error "Required file missing: $1"
}

check_dir() {
    [[ -d "$1" ]] || error "Required directory missing: $1"
}

cleanup_intermediate_files() {
    banner "Cleaning intermediate files"

    # Patterns to remove
    patterns=(
        "SPLITalnSeq*"
	"*mgaps"
        "*ntref"
	"*delta"
    )

    for pattern in "${patterns[@]}"; do
        found=false
        for f in $pattern; do
            if [[ -e "$f" ]]; then
                found=true
                rm -f "$f"
            fi
        done

        if [[ "$found" = true ]]; then
            info "Removed files matching: $pattern"
        else
            info "No files found for: $pattern"
        fi
    done

    info "Intermediate cleanup complete."
}

##############################
# Pre-flight Checks
##############################

banner "Checking required directories"

check_dir "$FASTA_DIR"
check_dir "$REF_DIR"
check_dir "$DATA_DIR"
check_dir "$ANNOT_DIR"
check_dir "$SCRIPTS_DIR"

banner "Checking required input files"

check_file "$MULTIFASTA"
check_file "$META_TABLE"

# Check segment-level reference fasta files, annotation files, and expected fasta outputs
for SEG in "${SEGMENTS[@]}"; do
    check_file "$REF_DIR/ref${SEG}.fa"
    check_file "$ANNOT_DIR/annot_table_${SEG}.pl"
done

info "All required directories and files are present."

# remove intermediate files
cleanup_intermediate_files

##############################
# Step 1: Split multifasta
##############################

banner "1. Splitting multifasta"
info "Creating segment-level FASTA files..."

perl "$SCRIPTS_DIR/reformatSequences.pl" "$MULTIFASTA"

info "Segment-level FASTA created."

##############################
# Step 2: Generate segment-level files
##############################

banner "2. Generating segment-level HaploCoV files"

for SEG in "${SEGMENTS[@]}"; do
    fasta_file="$DATA_DIR/${SEG}.segm.fa"
    check_file "$fasta_file"

    info "Processing segment: $SEG"

    perl "$SCRIPTS_DIR/addToTableNCBI.pl" \
        --metadata "$META_TABLE" \
        --ref "$REF_DIR/ref${SEG}.fa" \
        --seq "$fasta_file" \
	--nproc "${THREADS}" \
        --outfile "$DATA_DIR/${SEG}.HaploCoV"
done

info "All segment HaploCoV files generated."

##############################
# Step 3: Assign groups
##############################

banner "3. Assigning groups to HaploCoV files"

for FILE in "$DATA_DIR"/*.HaploCoV; do
    [[ -e "$FILE" ]] || error "No *.HaploCoV files found in data/"

    BASE=$(basename "$FILE")

    # Required auxiliary files
    check_file "${FILE}_V2_5_100"

    info "Assigning groups for: $BASE"

    perl "$SCRIPTS_DIR/assign.pl" \
        --dfile "${FILE}_V2_5_100" \
        --infile "$FILE" \
        --outfile "${FILE}.newD_v2_100"
done

info "All groups assigned."

##############################
# Step 4: Merge final table
##############################

banner "4. Merging all segment files"

perl "$SCRIPTS_DIR/mergeFilesFinal.pl" > "$OUTPUT_TABLE"

check_file "$OUTPUT_TABLE"

info "Merged table written to: $OUTPUT_TABLE"

##############################
# Step 5: Compute phenetic patterns
##############################

banner "5. Computing phenetic patterns"

perl "$SCRIPTS_DIR/DefToPhenetic.pl"

info "Phenetic patterns computed."

##############################
# Step 6: Annotation
##############################

banner "6. Annotating all segments"

for SEG in "${SEGMENTS[@]}"; do
    phenetic="$DATA_DIR/${SEG}.phenetic_V2.csv"

    info "Checking phenetic file for $SEG"
    check_file "$phenetic"

    info "Annotating: $SEG"

    perl "$SCRIPTS_DIR/annotate.pl" \
        --genome "$REF_DIR/ref${SEG}.fa" \
        --annot "$ANNOT_DIR/annot_table_${SEG}.pl" \
        --in "$phenetic" \
        --out "${SEG}.annot_V2.csv"
done

info "All segments annotated."

##############################
# Completed
##############################

banner "Pipeline finished successfully."

