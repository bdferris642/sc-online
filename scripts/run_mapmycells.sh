#!/usr/bin/env bash
set -euo pipefail

# defaults
clobber=0

usage() {
  echo "Usage: $0 -r REF.h5ad -q QUERY.h5ad -o OUTPUT_DIR -p PREFIX [-c]"
  echo "  -r, --ref           Path to reference .h5ad (required)"
  echo "  -q, --query         Path to query .h5ad (required)"
  echo "  -o, --output-path   Output directory (required)"
  echo "  -t, --hierarchy     Cell type hierarchy space-delimited (required)"
  echo "  -p, --prefix        Output name prefix (required)"
  echo "  -c, --clobber       Overwrite/recompute all intermediates (optional)"
  echo "  -h, --help          Show this help message and exit"
  exit 1
}

# Parse flags
hierarchy=()   # initialize as array

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref)
      ref="$2"
      shift
      ;;
    -q|--query)
      query="$2"
      shift
      ;;
    -o|--output-path)
      output_path="$2"
      shift
      ;;
    -p|--prefix)
      prefix="$2"
      shift
      ;;
    -t|--hierarchy)
      shift
      # collect all following non-flag arguments
      while [[ $# -gt 0 && "$1" != --* ]]; do
        hierarchy+=("$1")
        shift
      done
      continue   # IMPORTANT: skip the final shift below
      ;;
    -c|--clobber)
      clobber=1
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
  shift
done


# Required args
: "${ref:?Error: -r/--ref is required}"
: "${query:?Error: -q/--query is required}"
: "${output_path:?Error: -o/--output-path is required}"
: "${hierarchy:?Error: -t/--hierarchy is required}"
: "${prefix:?Error: -p/--prefix is required}"        # <-- enforce prefix

ref_dir="$(dirname "$ref")"

echo "Make output directory"
mkdir -p "$output_path"

chmod 777 "$ref_dir" || true
chmod 777 "$output_path" || true

echo "Reference file: $ref"
echo "Query file:     $query"
echo "Output dir:     $output_path"

# Build paths
ref_prefix="$(basename "$ref" .h5ad)"
ref_precomp_stats_path="${ref_dir}/${ref_prefix}_precomputed_stats.h5"
ref_marker_path="${ref_dir}/${ref_prefix}_reference_markers.h5"
query_marker_path="${output_path}/${prefix}_query_markers.json"   # <-- put in output dir

# hierarchy is a bash array: hierarchy=("cell_class" "cell_type")
json_hierarchy="[$(printf '"%s", ' "${hierarchy[@]}" | sed 's/, $//')]"
# -> '["cell_class", "cell_type"]'
echo "Using hierarchy: $json_hierarchy"

echo "== Precompute stats =="
if [[ ! -e "$ref_precomp_stats_path" || "$clobber" -eq 1 ]]; then
  echo "Computing stats for reference"
  python -m cell_type_mapper.cli.precompute_stats_scrattch \
    --h5ad_path "$ref" \
    --hierarchy "$json_hierarchy" \
    --n_processors 12 \
    --output_path "$ref_precomp_stats_path" \
    --clobber True \
    --normalization raw
else
  echo "Precomputed stats exist. Skipping."
fi

echo "== Reference markers =="
precompute_path_list="$(printf '["%s"]' "$ref_precomp_stats_path")"
if [[ ! -e "$ref_marker_path" || "$clobber" -eq 1 ]]; then
  echo "Computing reference markers"
  rm -f "${ref_dir}/reference_markers.h5"            # <-- force remove quietly
  python -m cell_type_mapper.cli.reference_markers \
    --precomputed_path_list "$precompute_path_list" \
    --n_valid 150 \
    --n_processors 36 \
    --output_dir "$ref_dir" \
    --clobber True

  # If tool salted the filename, prefer unsalted else pick newest salted file
  if [[ -e "${ref_dir}/reference_markers.h5" ]]; then
    mv "${ref_dir}/reference_markers.h5" "$ref_marker_path"
  else
    latest="$(ls -t "${ref_dir}"/reference_markers*.h5 2>/dev/null | head -n1 || true)"
    if [[ -n "$latest" ]]; then
      mv "$latest" "$ref_marker_path"
    else
      echo "Error: reference markers file not found." >&2
      exit 2
    fi
  fi
else
  echo "Reference markers exist. Skipping."
fi

echo "== Query markers =="
ref_marker_path_list="$(printf '["%s"]' "$ref_marker_path")"
if [[ ! -e "$query_marker_path" || "$clobber" -eq 1 ]]; then
  echo "Computing query markers"
  python -m cell_type_mapper.cli.query_markers \
    --output_path "$query_marker_path" \
    --reference_marker_path_list "$ref_marker_path_list" \
    --n_per_utility 100 \
    --n_processors 36
else
  echo "Query markers exist. Skipping."
fi

echo "== Label data =="
json_out="${output_path}/${prefix}_mapping_output.json"
csv_out="${output_path}/${prefix}_mapping_output.csv"
if [[ ! -e "$json_out" || ! -e "$csv_out" || "$clobber" -eq 1 ]]; then
  python -m cell_type_mapper.cli.from_specified_markers \
    --precomputed_stats.path "$ref_precomp_stats_path" \
    --query_markers.serialized_lookup "$query_marker_path" \
    --type_assignment.bootstrap_factor 0.5 \
    --type_assignment.bootstrap_iteration 1000 \
    --type_assignment.rng_seed 661123 \
    --type_assignment.n_processors 4 \
    --type_assignment.normalization raw \
    --query_path "$query" \
    --extended_result_path "$json_out" \
    --csv_result_path "$csv_out"
else
  echo "Mapping outputs exist. Skipping."
fi

echo "Done."
