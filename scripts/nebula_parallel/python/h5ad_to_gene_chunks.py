#!/usr/bin/env python
"""Stream a CSR count layer from an h5ad into gene-range chunks that R can load as dgCMatrix.

Replaces the "load full Seurat qs, then chunk-sobj-by-genes.R" step of run-nebula-parallel.sh
for objects whose total non-zeros exceed 2^31-1.

A CSR slice of cells x genes restricted to a contiguous gene range is, unchanged, the CSC
representation of genes x cells: p = per-cell pointers, i = local gene index, x = counts.
So each chunk is written as raw binary pieces per cell block:

    <raw-dir>/chunk_NNN/block_BBBBB.{n,i,x}.bin   (n: int32 per-cell entry counts,
                                                   i: int32 local gene index, x: float32 counts)
    <raw-dir>/chunk_NNN/genes.txt
    <raw-dir>/cells.txt, meta.tsv, meta_dtypes.json, chunk_plan.tsv, summary.json

write_chunk_qs.R then concatenates the blocks in order and builds one Seurat object per chunk.

Pass 1 computes per-gene nnz (for nnz-balanced chunk boundaries), per-cell total counts
(the default NEBULA offset, nCount_RNA) and checks counts are non-negative integers with
sorted column indices. Pass 2 writes the chunk pieces.
"""

import argparse
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from h5ad_io import (INT32_MAX, csr_info, ensure_env_python, plan_row_blocks, r_type_of,  # noqa: E402
                     read_csr_block, read_obs_derived, read_var_names, split_covariates)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pass1_block(args):
    """Per-gene nnz, per-cell sums/nnz and integrity checks for one row block."""
    h5ad, layer, r0, r1, n_var = args
    indptr, indices, data = read_csr_block(h5ad, layer, r0, r1)
    row_n = np.diff(indptr)
    rowid = np.repeat(np.arange(r1 - r0), row_n)
    gene_nnz = np.bincount(indices, minlength=n_var).astype(np.int64)
    row_sum = np.bincount(rowid, weights=data, minlength=r1 - r0)
    bad_values = int(np.count_nonzero((data < 0) | (data != np.round(data))))
    # indices must strictly increase within each row (row starts are allowed to drop)
    d = np.diff(indices)
    same_row = np.diff(rowid) == 0
    unsorted = int(np.count_nonzero((d <= 0) & same_row))
    out_of_range = int(np.count_nonzero((indices < 0) | (indices >= n_var)))
    return r0, gene_nnz, row_sum, row_n, bad_values, unsorted, out_of_range


def pass2_block(args):
    """Split one row block into per-chunk (n, i, x) binary pieces."""
    h5ad, layer, r0, r1, block_id, gene_starts, gene_ends, raw_dir = args
    indptr, indices, data = read_csr_block(h5ad, layer, r0, r1)
    n_rows = r1 - r0
    n_chunks = len(gene_starts)
    rowid = np.repeat(np.arange(n_rows, dtype=np.int64), np.diff(indptr))
    chunk_id = (np.searchsorted(gene_starts, indices, side="right") - 1).astype(np.int16)
    # stable sort keeps (row, gene) order inside each chunk -> valid CSC pieces
    order = np.argsort(chunk_id, kind="stable")
    counts_per_chunk = np.bincount(chunk_id, minlength=n_chunks)
    per_row = np.bincount(rowid * n_chunks + chunk_id, minlength=n_rows * n_chunks).reshape(n_rows, n_chunks)
    offsets = np.concatenate([[0], np.cumsum(counts_per_chunk)])
    for c in range(n_chunks):
        sel = order[offsets[c]:offsets[c + 1]]
        cdir = Path(raw_dir) / f"chunk_{c + 1:03d}"
        stem = cdir / f"block_{block_id:05d}"
        per_row[:, c].astype(np.int32).tofile(f"{stem}.n.bin")
        (indices[sel] - gene_starts[c]).astype(np.int32).tofile(f"{stem}.i.bin")
        data[sel].astype(np.float32).tofile(f"{stem}.x.bin")
    return block_id, counts_per_chunk


def nnz_balanced_ranges(gene_nnz, n_folds):
    """Contiguous gene ranges with roughly equal nnz. Returns (starts, ends) half-open, int64."""
    n_var = len(gene_nnz)
    if n_folds < 1 or n_folds > n_var:
        raise ValueError(f"--n-folds must be in [1, {n_var}], got {n_folds}")
    cum = np.cumsum(gene_nnz)
    targets = cum[-1] * np.arange(1, n_folds) / n_folds
    cuts = np.searchsorted(cum, targets, side="left") + 1
    # force strictly increasing cuts so that every chunk has >= 1 gene
    for k in range(len(cuts)):
        lo = (cuts[k - 1] + 1) if k > 0 else 1
        hi = n_var - (len(cuts) - k)
        cuts[k] = min(max(cuts[k], lo), hi)
    bounds = np.concatenate([[0], cuts, [n_var]]).astype(np.int64)
    return bounds[:-1], bounds[1:]


def export_meta(h5ad, raw_dir, id_col, covs_raw, offset_col, row_sum, derive_region_dummies):
    """Write meta.tsv / meta_dtypes.json for the columns NEBULA needs (id, covariates, nUMI_nebula)."""
    _, base_covs = split_covariates(covs_raw)
    # same columns as chunk-sobj-by-genes.R keeps: id, covariate parts, nUMI_nebula
    keep = list(dict.fromkeys([id_col, *base_covs]))
    from_obs = [c for c in keep if c != "nCount_RNA"]
    cols = list(dict.fromkeys(from_obs + ([] if offset_col == "nCount_RNA" else [offset_col])))
    obs = read_obs_derived(h5ad, cols, derive_region_dummies)
    if obs[id_col].isna().any():
        raise ValueError(f"id column '{id_col}' has {int(obs[id_col].isna().sum())} missing values")
    # 'nCount_RNA' always means the per-cell total of the count layer (what CreateSeuratObject
    # would compute on the full object); any obs column of that name is not used.
    if offset_col == "nCount_RNA":
        offset = row_sum
        log("offset nUMI_nebula = per-cell total of the count layer (nCount_RNA, full gene set)")
    else:
        offset = obs[offset_col].to_numpy(dtype=np.float64)
        log(f"offset nUMI_nebula = obs['{offset_col}']")
    if np.any(~np.isfinite(offset)) or np.any(offset <= 0):
        raise ValueError("offset must be finite and > 0 for every cell")
    meta = obs[from_obs].copy()
    # model.matrix silently drops rows with NA, which makes NEBULA's group_cell fail on every chunk
    n_na = meta.isna().sum()
    if n_na.any():
        raise ValueError("covariates with missing values (NEBULA cannot use them): "
                         + ", ".join(f"{c}={int(n)}" for c, n in n_na[n_na > 0].items()))
    if "nCount_RNA" in keep:
        meta["nCount_RNA"] = row_sum
    meta["nUMI_nebula"] = offset
    dtypes = {}
    for c in meta.columns:
        rt = r_type_of(meta[c])
        dtypes[c] = rt
        if rt == "logical":
            meta[c] = meta[c].map({True: "TRUE", False: "FALSE"})
        elif rt == "character":
            s = meta[c].astype(object)
            bad = s.dropna().astype(str).str.contains("\t|\n", regex=True)
            if bad.any():
                raise ValueError(f"column '{c}' contains tab/newline characters; cannot export as TSV")
            meta[c] = s
    meta.index.name = "cell"
    meta.to_csv(raw_dir / "meta.tsv", sep="\t", na_rep="NA")
    with open(raw_dir / "meta_dtypes.json", "w") as fh:
        json.dump(dtypes, fh, indent=2)
    for c, rt in dtypes.items():
        log(f"  meta column {c!r}: {rt}")
    return meta.index.to_numpy()


def main():
    ensure_env_python()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--layer", default="counts", help="CSR count layer ('X' for adata.X)")
    ap.add_argument("--raw-dir", required=True, help="output dir for raw chunk pieces")
    ap.add_argument("--n-folds", type=int, required=True)
    ap.add_argument("--id-col", required=True)
    ap.add_argument("--covs", default="")
    ap.add_argument("--offset-col", default="nCount_RNA",
                    help="'nCount_RNA' = computed per-cell total counts; otherwise an obs column")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--derive-region-dummies", action="store_true",
                    help="rebuild region_SN/region_VTA from region (SN->1,0; VTA->0,1; else 0,0)")
    ap.add_argument("--block-nnz", type=float, default=1e8, help="target stored entries per cell block")
    ap.add_argument("--min-free-gb", type=float, default=10.0)
    a = ap.parse_args()

    raw_dir = Path(a.raw_dir)
    if raw_dir.exists() and any(raw_dir.iterdir()):
        raise FileExistsError(f"{raw_dir} is not empty; remove it to re-chunk")
    raw_dir.mkdir(parents=True, exist_ok=True)

    n_obs, n_var, nnz, indptr = csr_info(a.h5ad, a.layer)
    log(f"{a.h5ad} layer={a.layer}: {n_obs:,} cells x {n_var:,} genes, nnz={nnz:,}")

    need_gb = nnz * 8 / 1e9 * 1.05 + n_obs * 4 * a.n_folds / 1e9
    free_gb = shutil.disk_usage(raw_dir).free / 1e9
    log(f"raw pieces need ~{need_gb:.1f} GB; free {free_gb:.1f} GB")
    if free_gb - need_gb < a.min_free_gb:
        raise OSError(f"Not enough disk: need ~{need_gb:.1f} GB + {a.min_free_gb} GB reserve, have {free_gb:.1f} GB")

    blocks = plan_row_blocks(indptr, int(a.block_nnz))
    log(f"{len(blocks)} cell blocks, {a.workers} workers")

    # ---- pass 1
    gene_nnz = np.zeros(n_var, dtype=np.int64)
    row_sum = np.zeros(n_obs, dtype=np.float64)
    row_nnz = np.zeros(n_obs, dtype=np.int64)
    bad = unsorted = oor = 0
    with ProcessPoolExecutor(a.workers) as ex:
        for r0, gn, rs, rn, b, u, o in ex.map(pass1_block, [(a.h5ad, a.layer, r0, r1, n_var) for r0, r1 in blocks]):
            gene_nnz += gn
            row_sum[r0:r0 + len(rs)] = rs
            row_nnz[r0:r0 + len(rn)] = rn
            bad, unsorted, oor = bad + b, unsorted + u, oor + o
    if bad or unsorted or oor:
        raise ValueError(f"Count layer failed checks: {bad} non-integer/negative values, "
                         f"{unsorted} unsorted/duplicate column indices, {oor} out-of-range indices")
    if gene_nnz.sum() != nnz:
        raise RuntimeError(f"pass 1 nnz {gene_nnz.sum()} != layer nnz {nnz}")
    log(f"pass 1 done: integer counts, sorted indices; {int((row_sum == 0).sum())} cells with zero counts")

    starts, ends = nnz_balanced_ranges(gene_nnz, a.n_folds)
    chunk_nnz = np.array([gene_nnz[s:e].sum() for s, e in zip(starts, ends)])
    if chunk_nnz.max() > INT32_MAX:
        raise ValueError(f"largest chunk has {chunk_nnz.max():,} nnz > 2^31-1; increase --n-folds "
                         f"(need at least {int(np.ceil(nnz / INT32_MAX * 1.2))})")
    log(f"{a.n_folds} gene chunks: nnz min {chunk_nnz.min():,} / median {int(np.median(chunk_nnz)):,} "
        f"/ max {chunk_nnz.max():,}; genes per chunk {int((ends - starts).min())}-{int((ends - starts).max())}")

    # ---- metadata (fails fast on missing columns before the long pass 2)
    cells = export_meta(a.h5ad, raw_dir, a.id_col, a.covs, a.offset_col, row_sum, a.derive_region_dummies)
    np.savetxt(raw_dir / "cells.txt", cells, fmt="%s")
    var_names = read_var_names(a.h5ad)
    if len(np.unique(var_names)) != n_var:
        raise ValueError("var index (gene names) is not unique")
    for c, (s, e) in enumerate(zip(starts, ends)):
        cdir = raw_dir / f"chunk_{c + 1:03d}"
        cdir.mkdir()
        np.savetxt(cdir / "genes.txt", var_names[s:e], fmt="%s")
    with open(raw_dir / "chunk_plan.tsv", "w") as fh:
        fh.write("chunk\tgene_start\tgene_end\tn_genes\tnnz\n")
        for c, (s, e, z) in enumerate(zip(starts, ends, chunk_nnz)):
            fh.write(f"{c + 1}\t{s + 1}\t{e}\t{e - s}\t{z}\n")

    # ---- pass 2
    written = np.zeros(a.n_folds, dtype=np.int64)
    jobs = [(a.h5ad, a.layer, r0, r1, b, starts, ends, str(raw_dir)) for b, (r0, r1) in enumerate(blocks)]
    with ProcessPoolExecutor(a.workers) as ex:
        for k, (_, cpc) in enumerate(ex.map(pass2_block, jobs), 1):
            written += cpc
            if k % max(1, len(jobs) // 10) == 0 or k == len(jobs):
                log(f"pass 2: {k}/{len(jobs)} blocks")
    if not np.array_equal(written, chunk_nnz):
        raise RuntimeError("pass 2 chunk nnz does not match pass 1 plan")

    summary = dict(h5ad=str(Path(a.h5ad).resolve()), layer=a.layer, n_cells=n_obs, n_genes=n_var, nnz=nnz,
                   n_folds=a.n_folds, n_blocks=len(blocks), id_col=a.id_col, covs=a.covs,
                   offset_col=a.offset_col, derive_region_dummies=a.derive_region_dummies, chunk_nnz=chunk_nnz.tolist())
    with open(raw_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    log("done")


if __name__ == "__main__":
    main()
