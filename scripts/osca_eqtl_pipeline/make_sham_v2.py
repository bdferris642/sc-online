#!/usr/bin/env python3
"""
make_sham_v2.py — Generate sham data and run v2 pipeline steps 2, 3, 5, 6, 7, 9.

Creates a self-contained test under sham_v2/ using:
  - First 60 samples of g1000_eur_hg38 as genotype data
  - 150 genes from gene_loc_v2.txt on chr1 as the gene universe
  - Simulated log-normal expression for 2 cell classes (Astrocyte, Microglia)
  - Simulated metadata (sex from FAM, random age/pmi/case_control)

Pipeline steps exercised:
  2  format OSCA inputs — checks ENSG IDs in Phenotype files and .opi probe column
  3  OSCA cis-eQTL
  5  two-stage FDR + RDS — checks v2 schema (ensg_id, gene_symbol; no Probe/Gene)
  6  common SNP-probe intersection (v2: SNP+ensg_id key)
  7  mashr
  9  validation — all 10 modules with ensg_id-keyed data
"""
from __future__ import annotations
import os, sys, shutil, subprocess
from pathlib import Path

# ── Self-activate into osca-venv ───────────────────────────────────────────────
SANDBOX   = Path(__file__).resolve().parents[1]
ENV_BIN   = SANDBOX / "micromamba_root" / "envs" / "osca-venv" / "bin"
if not ENV_BIN.is_dir():
    sys.exit("ERROR: osca-venv not found. Run setup.sh first.")
if sys.executable != str(ENV_BIN / "python"):
    os.execv(str(ENV_BIN / "python"), [str(ENV_BIN / "python")] + sys.argv)

import numpy as np
import pandas as pd
import pyreadr  # type: ignore

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = SANDBOX / "osca_eqtl_pipeline_v2"
SHAM_DIR   = SANDBOX / "sham_v2"
OSCA_BIN   = Path("/mnt/accessory/analysis/eqtl/osca/osca")
PLINK_DIR  = SANDBOX / "gwas" / "g1000_eur_hg38"
PLINK_PFX  = PLINK_DIR / "g1000_eur_hg38"
GENE_LOC   = SCRIPT_DIR / "gene_loc_v2.txt"
RESOURCES  = SANDBOX / "resources"

PLINK      = ENV_BIN / "plink2"
PYTHON     = ENV_BIN / "python"
RSCRIPT    = ENV_BIN / "Rscript"

# ── Parameters ─────────────────────────────────────────────────────────────────
N_SAMPLES    = 60
N_GENES      = 150        # genes from chr1 in gene_loc_v2
CELL_CLASSES = ["Astrocyte", "Microglia"]
N_THREADS    = 8
RNG_SEED     = 42

rng = np.random.default_rng(RNG_SEED)


def run(cmd: list, **kwargs) -> None:
    cmd = [str(c) for c in cmd]
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, **kwargs)


def run_sh(cmd: str, **kwargs) -> None:
    print(f"\n$ {cmd}")
    subprocess.run(cmd, shell=True, check=True, **kwargs)


# ── Directory setup ────────────────────────────────────────────────────────────
if SHAM_DIR.exists():
    shutil.rmtree(SHAM_DIR)

pb_dir       = SHAM_DIR / "pseudobulk"
osca_in_dir  = SHAM_DIR / "osca_inputs"
osca_out_dir = SHAM_DIR / "osca_outputs" / "eqtl_final_outs" / "sham"
val_dir      = SHAM_DIR / "validation"
for d in [pb_dir, osca_in_dir, osca_out_dir, val_dir]:
    d.mkdir(parents=True, exist_ok=True)


# ── 0. Load gene_loc_v2 — take first N_GENES genes on chr1 ───────────────────
print("\n=== 0. Loading gene_loc_v2 ===")
gene_loc = pd.read_csv(GENE_LOC, sep="\t", dtype={"ensg_id": str, "chr": str,
                                                    "gene_symbol": str, "strand": str})
gene_loc["TSS"] = gene_loc["TSS"].astype(int)
chr1_genes = gene_loc[gene_loc["chr"] == "1"].sort_values("TSS").head(N_GENES).reset_index(drop=True)
print(f"  Selected {len(chr1_genes)} genes from chr1")
print(f"  TSS range: {chr1_genes['TSS'].min():,} – {chr1_genes['TSS'].max():,}")


# ── 1. Subset PLINK bfiles ─────────────────────────────────────────────────────
print("\n=== 1. Subsetting PLINK bfiles ===")
plink_sub = SHAM_DIR / "geno_sham"
pca_pfx   = SHAM_DIR / "geno_sham_pca"

fam = pd.read_csv(str(PLINK_PFX) + ".fam", sep=r"\s+", header=None,
                  names=["FID","IID","PAT","MAT","SEX","PHENO"])
keep_file = SHAM_DIR / "keep_samples.txt"
fam.iloc[:N_SAMPLES][["FID","IID"]].to_csv(keep_file, sep="\t", index=False, header=False)

tss_min = chr1_genes["TSS"].min()
tss_max = chr1_genes["TSS"].max()
range_lo = max(1, tss_min - 1_000_000)
range_hi = tss_max + 1_000_000

# plink2 requires bim to be sorted by chr; sort first into pgen format then convert
plink_sorted = SHAM_DIR / "geno_sorted"
run([PLINK, "--bfile", PLINK_PFX, "--make-pgen", "--sort-vars", "--out", plink_sorted, "--silent"])
run([PLINK, "--pfile", plink_sorted,
     "--keep", keep_file,
     "--chr", "1", "--from-bp", range_lo, "--to-bp", range_hi,
     "--maf", "0.05", "--geno", "0.1",
     "--make-bed", "--out", plink_sub, "--silent"])

n_snps = sum(1 for _ in open(str(plink_sub) + ".bim"))
n_samp = sum(1 for _ in open(str(plink_sub) + ".fam"))
print(f"  Subset: {n_samp} samples × {n_snps:,} SNPs")

# PCA — plink2 writes a header line starting with #; strip it for plink1 compat
run([PLINK, "--bfile", plink_sub, "--pca", "5", "--out", pca_pfx, "--silent"])
eig_file = Path(str(pca_pfx) + ".eigenvec")
lines = eig_file.read_text().splitlines()
if lines and lines[0].startswith("#"):
    lines = lines[1:]
eig_file.write_text("\n".join(lines) + "\n")

# Copy bfiles + eigenvec into osca_in_dir (step 2 expects them there)
for ext in [".bed", ".bim", ".fam"]:
    shutil.copy(str(plink_sub) + ext, osca_in_dir / ("sham_vcf" + ext))
shutil.copy(str(pca_pfx) + ".eigenvec", osca_in_dir / "sham_vcf_pca.eigenvec")

# Step 2 R script hard-codes FID=0 in all phenotype/covariate files.
# OSCA matches individuals by FID+IID, so the PLINK FAM must also use FID=0.
sham_fam_path = osca_in_dir / "sham_vcf.fam"
sham_fam = pd.read_csv(sham_fam_path, sep=r"\s+", header=None,
                        names=["FID","IID","PAT","MAT","SEX","PHE"])
sham_fam["FID"] = 0
sham_fam.to_csv(sham_fam_path, sep="\t", header=False, index=False)

participants_file = SHAM_DIR / "participants.txt"
participants = fam.iloc[:N_SAMPLES]["IID"].tolist()
participants_file.write_text("\n".join(participants) + "\n")


# ── 2. Simulate expression + metadata per cell class ─────────────────────────
print("\n=== 2. Simulating pseudobulk expression + metadata ===")
fam_sub = fam.iloc[:N_SAMPLES].copy()

for cc in CELL_CLASSES:
    cc_safe = cc.replace(" ", "_")

    expr_mat = rng.lognormal(mean=2.5, sigma=0.8, size=(N_SAMPLES, N_GENES))
    expr_df = pd.DataFrame(expr_mat,
                           columns=chr1_genes["ensg_id"].tolist(),
                           index=fam_sub["IID"].tolist())
    # Write with unnamed index so R's read.csv() sees it as column "X" (matching step 1 output)
    expr_df.to_csv(pb_dir / f"{cc_safe}_expression_matrix_ds.csv")

    pd.DataFrame({"X": fam_sub["IID"].tolist(),
                  cc:  rng.uniform(0.05, 0.95, N_SAMPLES)
                  }).to_csv(pb_dir / f"{cc_safe}_composition_matrix_ds.csv", index=False)

    pd.DataFrame({
        "participant_id": fam_sub["IID"].tolist(),
        "sex":            fam_sub["SEX"].map({1:"M",2:"F"}).fillna("M").tolist(),
        "age":            rng.integers(50, 85, N_SAMPLES).astype(float),
        "pmi":            rng.uniform(5, 30, N_SAMPLES),
        "case_control":   rng.integers(0, 2, N_SAMPLES).astype(str),
        "study":          "sham",
        "brain_bank":     rng.choice(["BANK_A","BANK_B"], N_SAMPLES).tolist(),
        "dapi_nurr":      rng.uniform(0, 1, N_SAMPLES),
    }).to_csv(pb_dir / f"{cc_safe}_obs_metadata.csv", index=False)

print(f"  Written for: {CELL_CLASSES}")


# ── 3. Step 2: Format OSCA inputs ─────────────────────────────────────────────
print("\n=== 3. Step 2: Format OSCA inputs ===")
run([RSCRIPT, SCRIPT_DIR / "02-run-osca-formatting-scanpy.R",
     f"--expression-dir={pb_dir}",
     f"--output-dir={osca_in_dir}",
     "--vcf-slogan=sham_vcf",
     f"--participants={participants_file}",
     "--h5ad-cat-covars=sex case_control study brain_bank",
     "--h5ad-quant-covars=age pmi dapi_nurr",
     f"--gene-anot={GENE_LOC}"])

# Verify ENSG IDs in phenotype file and .opi probe column
print("\n  [verify step 2 outputs]")
for cc in CELL_CLASSES:
    cc_safe = cc.replace(" ", "_")
    pheno = osca_in_dir / f"Phenotype_{cc_safe}_osca.txt"
    opi   = osca_in_dir / f"Upprobe_{cc_safe}.opi"
    if not pheno.exists():
        print(f"  ERROR: {pheno.name} missing"); continue
    cols = pd.read_csv(pheno, sep="\t", nrows=0).columns.tolist()
    ensg_cols = [c for c in cols if c.startswith("ENSG")]
    bad_cols  = [c for c in cols if not c.startswith("ENSG") and c not in ("FID","IID")]
    opi_df    = pd.read_csv(opi, sep="\t", header=None, names=["chr","NAME","TSS","probe","strand"])
    ensg_probes = opi_df["probe"].str.startswith("ENSG").sum()
    print(f"  [{cc_safe}] Phenotype: {len(ensg_cols)} ENSG cols | non-ENSG gene cols: {bad_cols[:3] or 'none'}")
    print(f"  [{cc_safe}] .opi: {ensg_probes}/{len(opi_df)} probes are ENSG (expected {len(opi_df)})")


# ── 4. Step 3: OSCA cis-eQTL ──────────────────────────────────────────────────
print("\n=== 4. Step 3: OSCA cis-eQTL ===")
# Make osca available in PATH for 03-build-eqtl.sh
env = os.environ.copy()
env["PATH"] = str(OSCA_BIN.parent) + ":" + env.get("PATH", "")
env["LD_LIBRARY_PATH"] = "/mnt/accessory/anaconda3/lib" + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")

for cc in CELL_CLASSES:
    cc_safe = cc.replace(" ", "_")
    pheno   = osca_in_dir / f"Phenotype_{cc_safe}_osca.txt"
    befile  = osca_in_dir / f"befile_{cc_safe}"
    opi     = osca_in_dir / f"Upprobe_{cc_safe}.opi"
    cov1    = osca_in_dir / f"cov1_{cc_safe}.txt"
    cov2    = osca_in_dir / f"cov2_{cc_safe}_reduced.txt"
    out_tsv = osca_out_dir / f"eqtl_{cc_safe}"
    if not pheno.exists():
        print(f"  SKIP {cc_safe}: phenotype file missing"); continue
    print(f"  OSCA: {cc_safe} ...")
    run(["bash", SCRIPT_DIR / "03-build-eqtl.sh",
         pheno, befile, osca_in_dir / "sham_vcf",
         opi, cov1, cov2, N_THREADS, out_tsv], env=env)
    for f in osca_in_dir.glob(f"befile_{cc_safe}*"):
        f.rename(osca_out_dir / f.name)
    for f in osca_in_dir.glob(f"eqtl_{cc_safe}*log"):
        f.rename(osca_out_dir / f.name)
    for f in list(osca_in_dir.glob("tempeqtl_*")) + list(osca_out_dir.glob("tempeqtl_*")):
        f.unlink(missing_ok=True)


# ── 5. Step 5: FDR + RDS + plots ──────────────────────────────────────────────
print("\n=== 5. Step 5: FDR + RDS + plots ===")
# OSCA --query produces extensionless files (e.g. eqtl_Astrocyte, not eqtl_Astrocyte.tsv)
# Exclude .rds, .log, .bak.opi, .bod, .oii, .opi, .besd, .epi, .esi, and temp files
_EXCL_SUFFIXES = {".rds", ".log", ".opi", ".bod", ".oii", ".besd", ".epi", ".esi"}
osca_outputs = sorted(
    f for f in osca_out_dir.glob("eqtl_*")
    if f.is_file()
    and f.suffix not in _EXCL_SUFFIXES
    and ".bak" not in f.name
    and "tempeqtl_" not in f.name
    and "present_in_all" not in f.name
    and "sig_in_one" not in f.name
)
if not osca_outputs:
    print("  WARNING: no OSCA output files found — OSCA may have produced no output")
for osca_f in osca_outputs:
    print(f"  Processing {osca_f.name} ...")
    run([RSCRIPT, SCRIPT_DIR / "05-process-and-plot-osca-tsv.R", f"--path={osca_f}"])

print("\n  [verify step 5 RDS schema]")
errors = []
for rds_path in sorted(osca_out_dir.glob("eqtl_*.rds")):
    if "_sig" in rds_path.name or "present_in_all" in rds_path.name:
        continue
    r  = pyreadr.read_r(str(rds_path))
    df = list(r.values())[0]
    cols = list(df.columns)
    ok = "ensg_id" in cols and "gene_symbol" in cols and "Probe" not in cols and "Gene" not in cols
    print(f"  {rds_path.name}: ensg_id={'ensg_id' in cols} gene_symbol={'gene_symbol' in cols} "
          f"Probe={'Probe' in cols}(bad) Gene={'Gene' in cols}(bad) → {'OK' if ok else 'FAIL'}")
    if not ok:
        errors.append(rds_path.name)

# Remove OSCA raw output files (extensionless) after RDS are produced
for osca_f in osca_outputs:
    osca_f.unlink(missing_ok=True)


# ── 6. Step 6: Common SNP-probe intersection ───────────────────────────────────
print("\n=== 6. Step 6: Common SNP-probe intersection ===")
run([RSCRIPT, SCRIPT_DIR / "06-get-common-snp-probes-osca-rds.R",
     f"--base={osca_out_dir}"])


# ── 7. Step 7: mashr ──────────────────────────────────────────────────────────
print("\n=== 7. Step 7: mashr ===")
present_in_all = osca_out_dir / "eqtl_present_in_all.rds"
if present_in_all.exists():
    run([RSCRIPT, SCRIPT_DIR / "07-run-mashr.R",
         f"--path={present_in_all}",
         "--padj-thresh=0.05", "--num-random=1000", "--eps=1e-6"])
else:
    print("  SKIP: eqtl_present_in_all.rds not found")


# ── 8. Step 9: Validation (parallel per CC) ────────────────────────────────────
print("\n=== 8. Step 9: Validation ===")
procs = []
for cc in CELL_CLASSES:
    cc_safe = cc.replace(" ", "_")
    log_f = open(val_dir / f"{cc_safe}_step9.log", "w")
    p = subprocess.Popen(
        [str(PYTHON), str(SCRIPT_DIR / "09-eqtl-validation.py"),
         "--eqtl-dir",      str(osca_out_dir),
         "--gene-loc",      str(GENE_LOC),
         "--pb-output-dir", str(pb_dir),
         "--gtex-sn",       str(RESOURCES / "gtex_sn_signif_pairs.txt.gz"),
         "--atac-bed",      str(RESOURCES / "corces_2020_da_atac_peaks.bed.gz"),
         "--go-bp-gmt",     str(SANDBOX / "gene_sets/GO_Biological_Process_2025.gmt"),
         "--go-mf-gmt",     str(SANDBOX / "gene_sets/GO_Molecular_Function_2025.gmt"),
         "--out-dir",       str(val_dir),
         "--cell-class",    cc_safe],
        stdout=log_f, stderr=log_f
    )
    procs.append((cc_safe, p, log_f))
    print(f"  Launched {cc_safe} (pid {p.pid})")

val_failed = []
for cc_safe, p, log_f in procs:
    rc = p.wait(); log_f.close()
    print(f"  {cc_safe}: {'OK' if rc == 0 else f'FAILED (rc={rc})'}")
    if rc != 0:
        val_failed.append(cc_safe)
        lines = (val_dir / f"{cc_safe}_step9.log").read_text().splitlines()
        print("\n".join(lines[-30:]))


# ── Summary ────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("SHAM V2 PIPELINE SUMMARY")
print("="*60)
print(f"Output root: {SHAM_DIR}")
all_errors = errors + val_failed
for cc in CELL_CLASSES:
    cc_safe = cc.replace(" ", "_")
    cc_val  = val_dir / cc_safe
    n_plots = len(list(cc_val.glob("*.png"))) if cc_val.exists() else 0
    n_csvs  = len(list(cc_val.glob("*.csv"))) if cc_val.exists() else 0
    print(f"  {cc_safe}: {n_plots} plots, {n_csvs} CSVs in validation/")

if all_errors:
    print(f"\nFAILURES: {all_errors}")
    sys.exit(1)
else:
    print("\nAll steps passed.")
