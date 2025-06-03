#!/bin/bash

# packages:
# needs cell_type_mapper installed

# inputs:
# annotated h5ad path
# query h5ad path
# output path

# the cell_type_mapper.cli.precompute_stats_scrattch --hierarchy  is '["Class_label", "Subclass_label", "Group_label", "Cluster_label"]'

# Process the flags and arguments
while [[ "$#" -gt 0 ]]; do
	case $1 in
        	-r|--ref) # For flags for reference h5ad
	        	ref="$2"
			shift # Shift past the argument value
			;;
		-q|--query) # For flags for query h5ad
			query="$2"
			shift
			;;
		-o|--output_path) # 
			output_path="$2"
			shift
			;;
		-p|--prefix)
			prefix="$2"
			shift
			;;
		-c|--clobber)
			clobber=0
			;;
		-h|--help) # Help flag
			echo "Usage: $0 [-v|--verbose] [-f|--file filename] [-h|--help]"
			exit 0
			;;
		*) # Unknown option
			echo "Unknown option: $1"
			exit 1
			;;
	esac
	shift # Shift to the next argument
done

if [ -z "$ref" ]; then
	echo "Error: The -r or --ref flag option is required."
	exit 1
fi
if [ -z "$query" ]; then
	echo "Error: The -q or --query flag option is required."
	exit 1
fi
if [ -z "$output_path" ]; then
	echo "Error: The -o or --output_path flag option is required."
	exit 1
fi

ref_dir=$(dirname "${ref}")

echo "${ref_dir}"
echo "make output directory"
mkdir -p "${output_path}"

chmod 777 "${ref_dir}"
chmod 777 "${output_path}"

echo "Reference Directory: ${ref_dir}"
echo "Query Directory: ${query}"
echo "Output Directory: ${output_path}"

# Get precomputed stats
echo "get precomputed stats"
if [ ! -e "${ref_dir}/precomputed_stats.h5" ] || [ "${clobber}" -eq 1 ]; then
	echo "computing stats for reference"
	python -m cell_type_mapper.cli.precompute_stats_scrattch \
		--h5ad_path "${ref}" \
		--hierarchy '["cell_type_plot_condensed", "cell_type_plot"]' \
		--n_processors 12 \
		--output_path "${ref_dir}/precomputed_stats.h5" \
		--clobber True \
		--normalization raw
else
	echo "precomputed stats already exist. Skipping."
fi

# Get reference markers
echo "get reference markers"
precompute_path_list=$(printf '["%s"]' "${ref_dir}/precomputed_stats.h5")
if [ ! -e "${ref_dir}/reference_markers.h5" ] || [ "${clobber}" -eq 1 ]; then
	echo "computing reference markers"
	python -m cell_type_mapper.cli.reference_markers \
		--precomputed_path_list "$precompute_path_list" \
		--n_valid 30 \
		--n_processors 12 \
		--output_dir "${ref_dir}" \
		--clobber True
else
	echo "reference markers already exist. Skipping"
fi

# Get query markers
echo "get query markers"
ref_marker_path_list=$(printf '["%s"]' "${ref_dir}/reference_markers.h5")
if [ ! -e "${ref_dir}/query_markers.json" ] || [ "${clobber}" -eq 1 ]; then
	echo "computing query markers"
	python -m cell_type_mapper.cli.query_markers \
		--output_path "${ref_dir}/query_markers.json" \
		--reference_marker_path_list "$ref_marker_path_list" \
		--n_per_utility 100 \
		--n_processors 12
else
	echo "query markers already exist. Skipping"
fi

# Labeling data
echo "label data"
query_marker_path="${ref_dir}/query_markers.json"
if [ ! -e "${output_path}/mapping_output_higher_factor" ] || [ "${clobber}" -eq 1 ]; then
	python -m cell_type_mapper.cli.from_specified_markers \
		--precomputed_stats.path "${ref_dir}/precomputed_stats.h5" \
		--query_markers.serialized_lookup "${query_marker_path}" \
		--type_assignment.bootstrap_factor 0.5 \
		--type_assignment.bootstrap_iteration 1000 \
		--type_assignment.rng_seed 661123 \
		--type_assignment.n_processors 4 \
		--type_assignment.normalization raw \
		--query_path "$query" \
		--extended_result_path "${output_path}/${prefix}mapping_output.json" \
		--csv_result_path "${output_path}/${prefix}mapping_output.csv"
fi
