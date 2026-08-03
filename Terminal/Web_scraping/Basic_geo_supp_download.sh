#!/usr/bin/env bash

set -euo pipefail


ACCESSION="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"


########################################
# GEO series prefix
########################################

PREFIX=$(echo "$ACCESSION" | sed -E 's/^(GSE[0-9]{3}).*/\1nnn/')

USER_AGENT="Mozilla/5.0 (GEO downloader)"


########################################
# Validate file types
########################################

VALID_REGEX='\.(tar|tar\.gz|gz|csv|csv\.gz|rds|rds\.gz|h5ad|h5ad\.gz|txt|txt\.gz)$'


declare -A FILES=()


########################################
# Supplementary directory
########################################

SUPPL_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACCESSION}/suppl/"


echo
echo "Checking:"
echo "$SUPPL_URL"
echo


########################################
# Find files
########################################

while read -r filename; do

    if [[ "$filename" =~ $VALID_REGEX ]]; then

        FILES["$filename"]="${SUPPL_URL}${filename}"

    fi

done < <(

    curl \
        --fail \
        -s \
        -A "$USER_AGENT" \
        "$SUPPL_URL" \
    |
    grep -oE '[A-Za-z0-9._-]+\.(tar|tar.gz|gz|csv|csv.gz|rds|rds.gz|h5ad|h5ad.gz|txt|txt.gz)' \
    |
    sort -u

)



########################################
# Check results
########################################


echo "Found ${#FILES[@]} files"


if [[ ${#FILES[@]} -eq 0 ]]; then

    echo "ERROR: No supplementary files found"
    exit 1

fi


echo

for f in "${!FILES[@]}"; do
    echo "$f"
done

echo



########################################
# Download
########################################


for filename in "${!FILES[@]}"; do


    URL="${FILES[$filename]}"
    OUTPUT="${OUTDIR}/${filename}"


    echo "================================="
    echo "Downloading:"
    echo "$filename"
    echo "$URL"
    echo "================================="


    rm -f "${OUTPUT}.tmp"


    curl \
        --fail \
        -L \
        --http1.1 \
        -A "$USER_AGENT" \
        --retry 5 \
        --retry-delay 5 \
        -o "${OUTPUT}.tmp" \
        "$URL"



    ####################################
    # Check not HTML/XML error page
    ####################################


    FILETYPE=$(file "${OUTPUT}.tmp")


    if echo "$FILETYPE" | grep -qiE "HTML|XML|ASCII"; then

        echo
        echo "ERROR: Downloaded file is not data:"
        echo "$FILETYPE"
        echo

        head "${OUTPUT}.tmp"

        rm "${OUTPUT}.tmp"

        exit 1

    fi



    mv "${OUTPUT}.tmp" "$OUTPUT"



    ####################################
    # Validate gzip
    ####################################


    if [[ "$OUTPUT" == *.gz ]]; then

        echo "Checking gzip integrity..."

        gzip -t "$OUTPUT"

    fi



    echo "✓ Completed: $filename"
    echo


done



echo
echo "================================="
echo "DOWNLOAD COMPLETE"
echo "Saved to:"
echo "$OUTDIR"
echo "================================="
