#!/bin/bash
# Complete sequence download script
# Downloads all 7 required sequences from NCBI

set -e  # Exit on error

echo "========================================"
echo "SARS-CoV-2 Sequence Download Script"
echo "========================================"
echo ""

# Create directories
mkdir -p data/vaccine_references
mkdir -p data/sequences

# Counter
downloaded=0
failed=0

# Download function
download_sequence() {
    local url="$1"
    local output="$2"
    local name="$3"
    
    echo "Downloading $name..."
    
    if wget -q -O "$output" "$url"; then
        # Verify file is not empty and has content
        if [ -s "$output" ] && grep -q "^>" "$output"; then
            size=$(du -h "$output" | cut -f1)
            echo "  ✓ Downloaded $name ($size)"
            ((downloaded++))
        else
            echo "  ✗ Failed: Invalid file format"
            ((failed++))
        fi
    else
        echo "  ✗ Failed: Download error"
        ((failed++))
    fi
}

# Vaccine References
echo "--- Vaccine References ---"
download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=OR134577.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/vaccine_references/pfizer_bnt162b2_vector_OR134577.fasta" \
    "Pfizer (OR134577.1)"

download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=OR134578.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/vaccine_references/moderna_mrna1273_vector_OR134578.fasta" \
    "Moderna (OR134578.1)"

echo ""
echo "--- Natural Variants ---"

# Natural Sequences
download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=NC_045512.2&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta" \
    "Wuhan Reference (NC_045512.2)"

download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=MT020880.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/sequences/MT020880.1_Early_Wuhan.fasta" \
    "Early Wuhan (MT020880.1)"

download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=OM095706.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/sequences/OMX095706.1_Delta.fasta" \
    "Delta (OM095706.1)"

download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=OMX067679.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/sequences/OMX067679.1_Omicron_BA.1.fasta" \
    "Omicron BA.1 (OMX067679.1)"

download_sequence \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&log\$=seqview&db=nuccore&id=OMX067680.1&report=fasta&extrafeat=0&maxplex=1&hide-cdd=1" \
    "data/sequences/OMX067680.1_Omicron_BA.2.fasta" \
    "Omicron BA.2 (OMX067680.1)"

echo ""
echo "========================================"
echo "Download Summary"
echo "========================================"
echo "Downloaded: $downloaded/7"
echo "Failed: $failed/7"

if [ $downloaded -eq 7 ]; then
    echo ""
    echo "✅ All sequences downloaded successfully!"
    echo ""
    echo "Verifying downloads..."
    echo ""
    
    # Use seqkit if available, otherwise use grep
    if command -v seqkit &> /dev/null; then
        seqkit stats data/vaccine_references/*.fasta data/sequences/*.fasta
    else
        echo "seqkit not found. Basic verification:"
        for file in data/vaccine_references/*.fasta data/sequences/*.fasta; do
            if [ -f "$file" ]; then
                lines=$(wc -l < "$file")
                echo "$file: $lines lines"
            fi
        done
    fi
    
    echo ""
    echo "Ready to run analysis!"
    echo "Run: python3 independent_verification.py --data-dir data"
    exit 0
else
    echo ""
    echo "⚠ Some downloads failed. Please check your internet connection and try again."
    exit 1
fi
