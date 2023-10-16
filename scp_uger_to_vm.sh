#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <libnames.txt> <external_ip> <output_dir>"
        exit 1
fi

libnames_file="$1"
ip="$2"
output_dir="$3"

# If they don't exist, list the libraries in gs://macosko_data/libraries and save it to a file
if [ ! -e ~/macosko_data_libraries.txt]; then
    gsutil ls gs://macosko_data/libraries > macosko_data_libraries.txt
fi

# Loop through each library name in the specified file
while read -r name; do
    # Use grep to find the matching library path
    library_path=$(grep "$name" ~/macosko_data_libraries.txt)
    scp -r /Volumes/broad_macosko/data/multiome_demultiplex/vireo_output/"$name"/_GEX_vireo/ ferris@"$ip":~/"$output_dir"/$(basename "$library_path")
done < "$libnames_file"
