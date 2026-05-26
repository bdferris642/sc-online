#!/usr/bin/env python3
"""
run_gep_analysis.py — cNMF GEP analysis pipeline orchestrator.

Reads a manifest TSV and runs up to three analyses per row:
  regression   — variance partitioning / beta regression (run_gep_regression.R)
  intersection — Fisher's exact test gene-set overlap (get_cnmf_spectra_intersection.py)
  gsea         — GSEA (get_cnmf_spectra_gsea.py)

Usage:
    python cnmf-analysis/run_gep_analysis.py \\
        --manifest path/to/manifest.tsv \\
        --outdir   path/to/output/ \\
        [--steps   regression,intersection,gsea] \\
        [--jobs    N]

Manifest columns
----------------
Required:
  name              String identifier used as output subdirectory name (e.g. sox6-agtr1-k12)
  spectra_path      Path to *.spectra.k_*.dt_*.consensus.txt (GEPs × genes matrix)
  qs_path           Path to Seurat .qs object
  feature_loadings_path  Path to feature loadings file
  qc_geps           Comma-separated list of QC-passing GEP indices
  id_col            Column name for sample identifier
  covs              Comma-separated covariate column names

Optional (empty/NA → use defaults):
  ncores            Cores for regression (default: detectCores()-1)
  intersection_top_n     Max genes per GEP (default: 150)
  intersection_zscore    Z-score threshold (default: 3.0)
  intersection_min_size  Min gene-set size (default: 15)
  intersection_max_size  Max gene-set size (default: 250)
  gsea_jobs              Parallel GSEA jobs (default: n_geps)
  intersection_jobs      Parallel intersection jobs (default: n_geps)

Output layout
-------------
<outdir>/<name>/
  gep_regression/      ← run_gep_regression.R output
  intersection/        ← get_cnmf_spectra_intersection.py output
  gsea/                ← get_cnmf_spectra_gsea.py output
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR   = Path(__file__).parent
SANDBOX      = SCRIPT_DIR.parent
ENV_BIN      = SANDBOX / "micromamba_root" / "envs" / "cnmf-analysis-env" / "bin"
ANALYSIS_R   = ENV_BIN / "Rscript"
ANALYSIS_PY  = ENV_BIN / "python"
os.environ["LD_LIBRARY_PATH"] = (
    str(ENV_BIN.parent / "lib")
    + (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")
)
REGRESSION_R = SCRIPT_DIR / "run_gep_regression.R"
INTERSECTION_PY = SCRIPT_DIR / "get_cnmf_spectra_intersection.py"
GSEA_PY         = SCRIPT_DIR / "get_cnmf_spectra_gsea.py"

# Columns that go into the regression manifest (not pipeline-control columns)
REGRESSION_MANIFEST_COLS = [
    "qs_path", "feature_loadings_path", "qc_geps", "id_col", "covs",
]

# Optional columns that are pipeline-control only (never forwarded to sub-scripts)
PIPELINE_ONLY_COLS = {
    "name", "spectra_path",
    "ncores", "intersection_top_n", "intersection_zscore",
    "intersection_min_size", "intersection_max_size",
    "gsea_jobs", "intersection_jobs",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="cNMF GEP analysis pipeline: regression → intersection → GSEA."
    )
    parser.add_argument(
        "--manifest", type=Path, required=True,
        help="Path to manifest TSV"
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Base output directory"
    )
    parser.add_argument(
        "--steps", default="regression,intersection,gsea",
        help="Comma-separated steps to run (default: regression,intersection,gsea)"
    )
    parser.add_argument(
        "--jobs", type=int, default=1,
        help="Parallel manifest rows (default: 1, sequential)"
    )
    return parser.parse_args()


def load_manifest(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Manifest is empty: {path}")
    required = {"name", "spectra_path", "qs_path", "feature_loadings_path", "qc_geps", "id_col", "covs"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    return rows


def _val(row: dict, key: str) -> str | None:
    """Return non-empty, non-NA value or None."""
    v = row.get(key, "")
    if v is None or str(v).strip().lower() in ("", "na", "nan", "none"):
        return None
    return str(v).strip()


def write_regression_manifest(row: dict, tmp_path: Path):
    """Write a single-row TSV manifest for run_gep_regression.R."""
    cols = [c for c in REGRESSION_MANIFEST_COLS if c in row]
    with open(tmp_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerow({c: row.get(c, "") for c in cols})


def run_regression(row: dict, row_dir: Path):
    name = row["name"]
    out_subdir = row_dir / "gep_regression"
    out_subdir.mkdir(parents=True, exist_ok=True)

    if not ANALYSIS_R.exists():
        raise FileNotFoundError(
            f"cnmf-analysis-env Rscript not found: {ANALYSIS_R}\n"
            f"See cnmf-analysis/README.md for setup instructions."
        )

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", prefix=f"regression_manifest_{name}_",
        delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)

    try:
        write_regression_manifest(row, tmp_path)
        cmd = [
            str(ANALYSIS_R), "--no-save", "--no-restore", "--no-site-file", "--no-init-file",
            str(REGRESSION_R),
            "--input", str(tmp_path),
            "--output-dir", str(out_subdir),
        ]
        ncores = _val(row, "ncores")
        if ncores:
            cmd += ["--ncores", ncores]

        print(f"[{name}] Running regression: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"run_gep_regression.R failed for {name} (exit {result.returncode})")
    finally:
        tmp_path.unlink(missing_ok=True)

    print(f"[{name}] Regression complete → {out_subdir}")


def run_intersection(row: dict, row_dir: Path):
    name = row["name"]
    spectra_path = _val(row, "spectra_path")
    if not spectra_path:
        raise ValueError(f"[{name}] spectra_path is required for intersection step")

    out_subdir = row_dir / "intersection"

    cmd = [
        str(ANALYSIS_PY), str(INTERSECTION_PY),
        "--path", spectra_path,
        "--outdir", str(out_subdir),
    ]

    top_n    = _val(row, "intersection_top_n")
    zscore   = _val(row, "intersection_zscore")
    min_size = _val(row, "intersection_min_size")
    max_size = _val(row, "intersection_max_size")
    jobs     = _val(row, "intersection_jobs")

    if top_n:    cmd += ["--top-n",    top_n]
    if zscore:   cmd += ["--zscore",   zscore]
    if min_size: cmd += ["--min-size", min_size]
    if max_size: cmd += ["--max-size", max_size]
    if jobs:     cmd += ["--jobs",     jobs]

    print(f"[{name}] Running intersection: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"get_cnmf_spectra_intersection.py failed for {name} (exit {result.returncode})")

    print(f"[{name}] Intersection complete → {out_subdir}")


def run_gsea(row: dict, row_dir: Path):
    name = row["name"]
    spectra_path = _val(row, "spectra_path")
    if not spectra_path:
        raise ValueError(f"[{name}] spectra_path is required for gsea step")

    out_subdir = row_dir / "gsea"

    cmd = [
        str(ANALYSIS_PY), str(GSEA_PY),
        "--path", spectra_path,
        "--outdir", str(out_subdir),
    ]

    jobs = _val(row, "gsea_jobs")
    if jobs:
        cmd += ["--jobs", jobs]

    print(f"[{name}] Running GSEA: {' '.join(cmd)}")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"get_cnmf_spectra_gsea.py failed for {name} (exit {result.returncode})")

    print(f"[{name}] GSEA complete → {out_subdir}")


def process_row(row: dict, outdir: Path, steps: set[str]):
    name = row["name"]
    row_dir = outdir / name
    row_dir.mkdir(parents=True, exist_ok=True)

    if "regression" in steps:
        run_regression(row, row_dir)
    if "intersection" in steps:
        run_intersection(row, row_dir)
    if "gsea" in steps:
        run_gsea(row, row_dir)


def main():
    args = parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {args.manifest}")

    steps = {s.strip().lower() for s in args.steps.split(",")}
    valid_steps = {"regression", "intersection", "gsea"}
    unknown = steps - valid_steps
    if unknown:
        raise ValueError(f"Unknown steps: {sorted(unknown)}. Valid: {sorted(valid_steps)}")

    rows = load_manifest(args.manifest)
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Manifest: {args.manifest} ({len(rows)} row(s))")
    print(f"Steps:    {', '.join(sorted(steps))}")
    print(f"Output:   {args.outdir}")
    print(f"Jobs:     {args.jobs}\n")

    if args.jobs == 1:
        errors = []
        for row in rows:
            try:
                process_row(row, args.outdir, steps)
            except Exception as exc:
                errors.append((row["name"], exc))
                print(f"ERROR [{row['name']}]: {exc}", file=sys.stderr)
    else:
        errors = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(process_row, row, args.outdir, steps): row["name"]
                for row in rows
            }
            for fut in as_completed(futures):
                name = futures[fut]
                try:
                    fut.result()
                except Exception as exc:
                    errors.append((name, exc))
                    print(f"ERROR [{name}]: {exc}", file=sys.stderr)

    if errors:
        print(f"\n{len(errors)} row(s) failed:", file=sys.stderr)
        for name, exc in errors:
            print(f"  {name}: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"\nAll rows complete. Results in: {args.outdir}")


if __name__ == "__main__":
    main()
