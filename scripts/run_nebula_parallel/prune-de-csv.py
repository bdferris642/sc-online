import argparse
import pandas as pd

COLS = [
    "gene",
    "logFC_case_controlpd",
    "padj_case_controlpd",
    "rank_purity",
    "signed_neg_log10_p_case_controlpd",
    "p_case_controlpd",
    "signed_neg_log10_padj_case_controlpd",
    "se_case_controlpd",
    "logFC_(Intercept)"
]
GENE_COL = "gene"
GENE_ID_COL="gene_id"
ensg_to_symbol_map = pd.read_csv("/home/ferris/sc-online/misc/ensg_to_symbol.csv")

parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True, help='Path to the CSV file')
args = parser.parse_args() 

df = pd.read_csv(args.path, index_col=0)
df = df[COLS]

if all(df[GENE_COL].str.startswith("ENSG")):
    df[GENE_ID_COL] = df[GENE_COL]
    df = pd.merge(df, ensg_to_symbol_map, left_on=GENE_ID_COL, right_on="ensembl_gene_id")
    df[GENE_COL] = df["hgnc_symbol"]
    df = df[COLS]

df["significant"] = df["padj_case_controlpd"] < 0.05

# order by signed_neg_log10_padj_case_controlpd desc
df = df.sort_values("signed_neg_log10_padj_case_controlpd", ascending=False)

df.to_csv(args.path.replace(".csv", ".pruned.csv"))