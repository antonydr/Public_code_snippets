#!/usr/bin/env bash

set -euo pipefail

ACCESSION="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"

PREFIX="${ACCESSION:0:6}nnn"

USER_AGENT="GEO-Downloader/1.0"

VALID_REGEX='\.(tar|tar\.gz|gz|csv|csv\.gz|rds|rds\.gz|h5ad|h5ad\.gz|txt|txt\.gz)$'

declare -A FILES=()


########################################
# Normalize FTP URLs to HTTPS
########################################

normalize_url() {

    local url="$1"

    if [[ "$url" == ftp://ftp.ncbi.nlm.nih.gov* ]]; then
        url="${url/ftp:\/\//https:\/\/}"
    fi

    echo "$url"
}


########################################
# Add unique downloadable file
########################################

add_file() {

    local url="$1"

    url=$(normalize_url "$url")

    local filename
    filename=$(basename "$url")

    [[ -z "$filename" ]] && return


    if [[ "$filename" =~ $VALID_REGEX ]]; then

        if [[ "$url" =~ (biorxiv|doi.org|pubmed) ]]; then
            return
        fi

        FILES["$filename"]="$url"

    fi
}


########################################
# 1. GEO accession page extraction
########################################

echo "Extracting GEO links..."

ACC_URL="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${ACCESSION}"


while read -r link; do

    add_file "$link"

done < <(

    curl \
        -s \
        -H "User-Agent: ${USER_AGENT}" \
        "$ACC_URL" \
    | grep -oE 'href="[^"]+"' \
    | sed 's/href="//;s/"//'

)


########################################
# 2. GEO supplementary directory
########################################

SUPPL_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACCESSION}/suppl/"

echo "Checking supplement directory:"
echo "$SUPPL_URL"


while read -r filename; do

    add_file "${SUPPL_URL}${filename}"

done < <(

    curl \
        -s \
        -H "User-Agent: ${USER_AGENT}" \
        "$SUPPL_URL" \
    | grep -oE '[A-Za-z0-9._-]+\.(tar|tar.gz|gz|csv|csv.gz|rds|rds.gz|h5ad|h5ad.gz|txt|txt.gz)' \
    | sort -u

)


########################################
# 3. GEO text listing fallback
########################################

TEXT_URL="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${ACCESSION}&targ=self&view=full&form=text"


while read -r link; do

    add_file "$link"

done < <(

    curl \
        -s \
        -H "User-Agent: ${USER_AGENT}" \
        "$TEXT_URL" \
    | grep "$ACCESSION" \
    | grep -oE '[^[:space:]]+\.(tar|tar.gz|gz|csv|csv.gz|rds|rds.gz|h5ad|h5ad.gz|txt|txt.gz)' \
    || true

)


########################################
# Summary
########################################

echo
echo "Found ${#FILES[@]} unique downloadable files"
echo


if [[ ${#FILES[@]} -eq 0 ]]; then

    echo "ERROR: No downloadable files found"
    exit 1

fi


INDEX=1

for filename in "${!FILES[@]}"; do

    echo "  [$INDEX] $filename"
    INDEX=$((INDEX+1))

done


echo


########################################
# Download files
########################################

TOTAL=${#FILES[@]}
COUNT=0


for filename in "${!FILES[@]}"; do

    COUNT=$((COUNT+1))

    URL="${FILES[$filename]}"
    OUTPUT="${OUTDIR}/${filename}"


    echo "========================================"
    echo "[$COUNT/$TOTAL] Downloading:"
    echo "$filename"
    echo "$URL"
    echo "========================================"


    if [[ -f "$OUTPUT" ]]; then
        echo "Existing file found, resuming if incomplete..."
    fi


    curl \
        -L \
        --http1.1 \
        -H "User-Agent: ${USER_AGENT}" \
        -H "Accept: */*" \
        --retry 5 \
        --retry-delay 5 \
        --continue-at - \
        -o "$OUTPUT" \
        "$URL"


    echo
    echo "✓ Finished: $filename"


    ####################################
    # SHA256 checksum
    ####################################

    if command -v sha256sum >/dev/null 2>&1; then

        SHA=$(sha256sum "$OUTPUT" | awk '{print $1}')

    else

        SHA=$(shasum -a 256 "$OUTPUT" | awk '{print $1}')

    fi


    echo "SHA256: $SHA"
    echo

done


echo "========================================"
echo "DOWNLOAD COMPLETE"
echo "Files saved to:"
echo "$OUTDIR"
echo "========================================"
