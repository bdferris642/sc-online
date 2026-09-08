#!/usr/bin/env bash
# 08-run-smr.sh — Run SMR + HEIDI for each cell class against a GWAS trait.
#
# SMR (Summary-based Mendelian Randomization) tests whether a GWAS signal is
# mediated by a cis-eQTL effect on gene expression. For each top cis-eQTL SNP
# in a gene's locus, it asks: is the GWAS association at that SNP explained by
# the gene's expression? The HEIDI test then rejects loci where the pattern is
# better explained by two independent causal SNPs in LD rather than one shared
# causal variant (pleiotropy vs. linkage).
#
# Pipeline per cell class:
#   0. Build rs-only bfile from the LD reference panel (once; written to out-dir)
#   1. Locate or create SMR qfile TSV (priority: out-dir > eqtl-dir > re-export RDS)
#   2. Create BESD (SMR binary eQTL format) from qfile (if not already present)
#   3. Run SMR + HEIDI against GWAS summary statistics
#
# Required arguments:
#   --eqtl-dir DIR       Directory containing eqtl_{cc}.rds files (step 5 output)
#   --bfile PREFIX       PLINK bfile prefix for LD reference panel (e.g. sn_vta_annot)
#                        Must contain only rsIDs; a filtered copy is auto-created if needed.
#   --gwas FILE          GWAS summary stats in SMR .ma format (SNP A1 A2 freq b se p N)
#   --out-dir DIR        Output directory for SMR results
#
# Optional arguments:
#   --smr-bin PATH       Path to SMR binary              [auto-detected from PATH]
#   --rscript PATH       Path to Rscript                  [auto-detected from PATH]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"

# ── Locate env ────────────────────────────────────────────────────────────────
_ENV1="${SANDBOX}/micromamba_root/envs/osca-venv"
_ENV2="/home/ferris/sc-online/scripts/micromamba_root/envs/osca-venv"
if   [ -d "${_ENV1}/bin" ]; then ENV_BIN="${_ENV1}/bin"
elif [ -d "${_ENV2}/bin" ]; then ENV_BIN="${_ENV2}/bin"
else echo "ERROR: osca-venv not found. Run setup.sh first."; exit 1; fi

# ── Defaults ──────────────────────────────────────────────────────────────────
EQTL_DIR=""
BFILE=""
GWAS=""
OUT_DIR=""
SMR_BIN="${ENV_BIN}/smr"
RSCRIPT="${ENV_BIN}/Rscript"

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --eqtl-dir)  EQTL_DIR="$2";  shift ;;
        --bfile)     BFILE="$2";     shift ;;
        --gwas)      GWAS="$2";      shift ;;
        --out-dir)   OUT_DIR="$2";   shift ;;
        --smr-bin)   SMR_BIN="$2";   shift ;;
        --rscript)   RSCRIPT="$2";   shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

: "${EQTL_DIR:?--eqtl-dir is required}"
: "${BFILE:?--bfile is required}"
: "${GWAS:?--gwas is required}"
: "${OUT_DIR:?--out-dir is required}"

# Fall back to system smr if not in env bin
if [ ! -f "${SMR_BIN}" ]; then
    SMR_BIN="$(command -v smr 2>/dev/null || true)"
    if [ -z "${SMR_BIN}" ]; then
        # Try the known location
        SMR_BIN="/mnt/accessory/analysis/eqtl/smr-1.3.1-linux-x86_64/smr"
    fi
fi
if [ ! -f "${SMR_BIN}" ]; then
    echo "ERROR: smr binary not found. Pass --smr-bin or add to PATH."
    exit 1
fi
echo "Using SMR: ${SMR_BIN}"

if [ ! -f "${GWAS}" ]; then echo "ERROR: GWAS file not found: ${GWAS}"; exit 1; fi
if [ ! -f "${BFILE}.bed" ]; then echo "ERROR: bfile not found: ${BFILE}.bed"; exit 1; fi

mkdir -p "${OUT_DIR}"

# ── Step 0: create rs-only bfile (SMR requires rsIDs; done once) ──────────────
# SMR will error on any SNP whose ID does not start with "rs".
# We extract only rsID SNPs from the LD reference bfile and write the filtered
# copy into out-dir (which is always writable, unlike the source bfile location).
BFILE_RS="${OUT_DIR}/$(basename "${BFILE}")_rs_only"
if [ ! -f "${BFILE_RS}.bed" ]; then
    echo "==> Step 0: building rs-only bfile ..."
    echo "    Source  : ${BFILE}.{bed,bim,fam}"
    echo "    Output  : ${BFILE_RS}.{bed,bim,fam}"
    # Count total and rs SNPs for user visibility
    total_snps=$(wc -l < "${BFILE}.bim")
    rs_snps=$(awk '$2 ~ /^rs/' "${BFILE}.bim" | wc -l)
    echo "    SNPs in source bfile : ${total_snps}  (rs-only: ${rs_snps})"
    awk '$2 ~ /^rs/{print $2}' "${BFILE}.bim" > "${OUT_DIR}/rs_snplist.txt"
    "${ENV_BIN}/plink2" \
        --bfile   "${BFILE}" \
        --extract "${OUT_DIR}/rs_snplist.txt" \
        --make-bed \
        --out     "${BFILE_RS}" \
        --threads 50
    echo "    rs-only bfile written to ${BFILE_RS}"
else
    echo "==> Step 0: rs-only bfile already exists — ${BFILE_RS}"
fi

# ── Discover cell classes from per-cell-class RDS files ───────────────────────
mapfile -t RDS_FILES < <(ls "${EQTL_DIR}"/eqtl_*.rds 2>/dev/null | grep -v "_sig\.rds" | grep -v "present_in_all")
if [ ${#RDS_FILES[@]} -eq 0 ]; then
    echo "ERROR: No eqtl_*.rds files found in ${EQTL_DIR}"
    exit 1
fi
echo "Cell class RDS files found:"
printf '  %s\n' "${RDS_FILES[@]}"

# ── Process each cell class ───────────────────────────────────────────────────
for RDS in "${RDS_FILES[@]}"; do
    CC=$(basename "${RDS}" .rds | sed 's/^eqtl_//')
    QFILE="${OUT_DIR}/eqtl_${CC}.tsv"
    PREFIX="${OUT_DIR}/eqtl_${CC}"
    echo ""
    echo "══════════════════════════════════════════"
    echo " Cell class: ${CC}"
    echo "══════════════════════════════════════════"

    # ── Step 1: build rs-only SMR qfile TSV ───────────────────────────────────
    # SMR qfile format: SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b se p
    # (note lowercase 'se'; OSCA outputs 'SE' — normalised on read)
    # Always regenerated so rs-only filtering is guaranteed; BESD creation (step 2)
    # is the idempotent checkpoint that avoids redundant work on re-runs.
    # Source priority: OSCA TSV in eqtl-dir (fast) → RDS (slower, full reload).
    EQTL_TSV="${EQTL_DIR}/eqtl_${CC}.tsv"
    if [ -f "${EQTL_TSV}" ]; then
        echo "[${CC}] Step 1: filtering OSCA TSV → rs-only qfile ..."
        "${RSCRIPT}" --no-save --no-restore -e "
            suppressMessages(library(dplyr))
            tsv <- read.table('${EQTL_TSV}', header=TRUE, sep='\t', stringsAsFactors=FALSE)
            qfile <- tsv %>%
                filter(grepl('^rs', SNP)) %>%
                select(SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr, Probe_bp,
                       Gene, Orientation, b, se=SE, p)
            write.table(qfile, '${QFILE}', sep='\t', quote=FALSE, row.names=FALSE)
            cat(sprintf('  Written %d rows to ${QFILE}\n', nrow(qfile)))
        "
    else
        echo "[${CC}] Step 1: no TSV in eqtl-dir — exporting from RDS ..."
        "${RSCRIPT}" --no-save --no-restore -e "
            suppressMessages(library(dplyr))
            rds <- readRDS('${RDS}')
            qfile <- rds %>%
                filter(grepl('^rs', SNP)) %>%
                select(SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr, Probe_bp,
                       Gene, Orientation, b, se=SE, p)
            write.table(qfile, '${QFILE}', sep='\t', quote=FALSE, row.names=FALSE)
            cat(sprintf('  Written %d rows to ${QFILE}\n', nrow(qfile)))
        "
    fi

    # ── Step 2: create BESD (SMR binary eQTL format) ──────────────────────────
    if [ ! -f "${PREFIX}.besd" ]; then
        echo "[${CC}] Creating BESD file ..."
        "${SMR_BIN}" --qfile "${QFILE}" --make-besd --out "${PREFIX}"
    else
        echo "[${CC}] BESD already exists: ${PREFIX}.besd"
    fi

    # ── Step 3: run SMR + HEIDI ───────────────────────────────────────────────
    echo "[${CC}] Running SMR + HEIDI ..."
    "${SMR_BIN}" \
        --bfile     "${BFILE_RS}" \
        --gwas-summary "${GWAS}" \
        --beqtl-summary "${PREFIX}" \
        --out       "${PREFIX}_smr"

    echo "[${CC}] Done → ${PREFIX}_smr.smr"
done

echo ""
echo "All cell classes complete."
