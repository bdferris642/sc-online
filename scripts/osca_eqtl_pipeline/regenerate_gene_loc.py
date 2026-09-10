#!/usr/bin/env python
# regenerate_gene_loc.py — Rebuild gene_loc.txt from MyGene.info using Ensembl IDs.
#
# Queries MyGene.info in batches for all Ensembl IDs in the expression data,
# retrieves: Entrez ID (probe), chromosome, TSS, HGNC symbol, strand.
# Filters to chromosomes 1-22 and writes OSCA-compatible gene_loc.txt.
#
# Usage:
#   python regenerate_gene_loc.py \
#     --expr-csv  /path/to/any_expression_matrix_ds.csv \
#     --out       /path/to/gene_loc_new.txt

import argparse
import time
import pandas as pd
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--expr-csv", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--batch-size", type=int, default=1000)
args = parser.parse_args()

# ── 1. Get Ensembl IDs from expression CSV header ─────────────────────────────
print("Reading Ensembl IDs from expression CSV ...")
header = pd.read_csv(args.expr_csv, nrows=0)
ensg_ids = [c for c in header.columns if c.startswith("ENSG")]
print(f"  {len(ensg_ids)} Ensembl IDs found")

# ── 2. Query MyGene.info in batches ───────────────────────────────────────────
FIELDS = "entrezgene,symbol,genomic_pos,genomic_pos_hg19"
URL    = "https://mygene.info/v3/gene"

results = []
for i in range(0, len(ensg_ids), args.batch_size):
    batch = ensg_ids[i:i + args.batch_size]
    print(f"  Querying batch {i//args.batch_size + 1} / {-(-len(ensg_ids)//args.batch_size)} ...")
    for attempt in range(3):
        try:
            r = requests.post(URL,
                data={"ids": ",".join(batch), "fields": FIELDS,
                      "species": "human", "dotfield": "false"},
                timeout=60)
            r.raise_for_status()
            results.extend(r.json())
            break
        except Exception as e:
            if attempt == 2:
                raise
            print(f"    Retrying after error: {e}")
            time.sleep(5)

print(f"  Retrieved {len(results)} records")

# ── 3. Parse into gene_loc rows ───────────────────────────────────────────────
AUTOSOMES = {str(c) for c in range(1, 23)}

rows = []
for rec in results:
    if rec.get("notfound"):
        continue

    entrez = rec.get("entrezgene")
    symbol = rec.get("symbol", "")
    gpos   = rec.get("genomic_pos")

    if not entrez or not symbol or not gpos:
        continue

    # genomic_pos can be a list (multiple locations) or a single dict
    locs = gpos if isinstance(gpos, list) else [gpos]

    for loc in locs:
        chrom = str(loc.get("chr", ""))
        if chrom not in AUTOSOMES:
            continue
        start  = loc.get("start")
        end    = loc.get("end")
        strand = loc.get("strand")   # 1 or -1
        if start is None or strand is None:
            continue
        # TSS = start for + strand, end for - strand
        tss = start if strand == 1 else end
        rows.append({
            "probe":  int(entrez),
            "chr":    int(chrom),
            "TSS":    int(tss),
            "NAME":   symbol,
            "strand": "+" if strand == 1 else "-"
        })

df = pd.DataFrame(rows)
if df.empty:
    raise RuntimeError("No rows parsed — check MyGene.info response format")

# One row per symbol: pick the canonical location (most 5' TSS)
df = (df.sort_values(["NAME", "TSS"])
        .drop_duplicates(subset=["NAME"], keep="first")
        .drop_duplicates(subset=["probe"], keep="first")
        .sort_values(["chr", "TSS"])
        .reset_index(drop=True))

print(f"Final gene_loc: {len(df)} genes across chromosomes 1-22")

# ── 4. Write output ───────────────────────────────────────────────────────────
import os
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
df.to_csv(args.out, sep="\t", index=False)
print(f"Written to: {args.out}")

# ── 5. Coverage report ────────────────────────────────────────────────────────
queried = set(ensg_ids)
found   = {r["query"] for r in results if not r.get("notfound")}
print(f"\nCoverage:")
print(f"  Queried:         {len(queried)}")
print(f"  Found by MGI:    {len(found)} ({100*len(found)/len(queried):.1f}%)")
print(f"  In gene_loc:     {len(df)}")
