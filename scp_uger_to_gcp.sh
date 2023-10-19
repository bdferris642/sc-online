#!/bin/zsh

# Ensure we have the right number of arguments
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <path_to_input_names_file> <ip> <output_dir_on_gcloud>"
    exit 1
fi

INPUT_NAMES_FILE="$1"
IP="$2"
OUTPUT_DIR="$3"

rm missing_from_vireo.txt 

# For every name in the INPUT_NAMES
while IFS= read -r name; do
    echo ::::::::::::::::::::::::::::::::::::::::::::::::
    # Use SSH to make the directory in ferris@<IP>:OUTPUT_DIR/name
    ssh ferris@"$IP" "mkdir -p $OUTPUT_DIR/$name" < /dev/null

    # Define a LOCAL_DIR by grep'ing the name from the INPUT_NAMES in the file `vireo_output.txt`
    LOCAL_DIR=/Volumes/broad_macosko/data/multiome_demultiplex/vireo_output/"$(grep $name vireo_output.txt)"
    echo "$LOCAL_DIR"

    # If it doesn't exist or if there are multiple matches
    if [[ -z "$LOCAL_DIR" ]]; then
        echo "$name" >> missing_from_vireo.txt
        continue
    elif [[ $(echo "$LOCAL_DIR" | wc -l) -gt 1 ]]; then
        # If there's more than one LOCAL_DIR, take the longest name that does not contain the prefix "HOLD"
        LOCAL_DIR=$(echo "$LOCAL_DIR" | grep -v "HOLD" | awk '{print length, $0}' | sort -n | tail -n 1 | cut -d' ' -f2-)
    fi

    # Check for _GEX_vireo contents
    if [[ -n "$(ls $LOCAL_DIR/_GEX_vireo 2>/dev/null)" ]]; then
        LOCAL_SUBDIR="$LOCAL_DIR/_GEX_vireo"
    # Check for yymmdd_GEX_vireo contents
    elif [[ -n "$(ls $LOCAL_DIR/[0-9][0-9][0-9][0-9][0-9][0-9]_GEX_vireo 2>/dev/null)" ]]; then
        LOCAL_SUBDIR="$LOCAL_DIR/$(ls $LOCAL_DIR | grep '[0-9][0-9][0-9][0-9][0-9][0-9]_GEX_vireo')"
    else
        echo "$name" >> missing_from_vireo.txt
        continue
    fi

    echo "$LOCAL_SUBDIR"

    # Use SCP to copy the CONTENTS of LOCAL_SUBDIR to ferris@<IP>:OUTPUT_DIR/name
    scp $LOCAL_SUBDIR/* ferris@"$IP":"$OUTPUT_DIR/$name/" < /dev/null

    # Unpack the .gz files on the remote server
    ssh ferris@"$IP" "gzip -d $OUTPUT_DIR/$name/*.gz" < /dev/null

    echo ::::::::::::::::::::::::::::::::::::::::::::::::

done < "$INPUT_NAMES_FILE"

