#!/usr/bin/env python3
"""
simulate_data.py
Simulate a cell-level snRNA-seq CSV based on parameters from
data-models/sn-vta-neuron-obs-info.json.

Design:
  - 40 donors (20 PD, 20 CTR), most with 2 tissue-region samples (SN + VTA)
    so that donor_id is a meaningful random effect.
  - Cell type proportions derived from the 33-type classification.
  - In PD donors: Dopaminergic neurons depleted ~45%.
  - In PD donors: Da_SOX6_EYA4 depleted ~75% (strongest depletion).
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# ── Cell type roster (matches sn-vta-neuron-obs-info.json) ────────────────────
CELL_TYPES = [
    "Da_CALB1_CALCR", "Da_CALB1_CRYM", "Da_CALB1_TMEM200A_GAD2",
    "Da_CALB1_TMEM200A_GEM", "Da_CALB1_TMEM200A_TRHR",
    "Da_SOX6_AGTR1", "Da_SOX6_AGTR1_EYA4_NXPH1", "Da_SOX6_EYA4",
    "Da_SOX6_GFRA2_PART1", "Da_SOX6_NXPH1",
    "Exc_5HT", "Exc_CYP2J2_CRH", "Exc_HDC", "Exc_LHX9_FBN2",
    "Exc_LMX1A_PITX2", "Exc_NEUROD1", "Exc_PAX5_F13A1", "Exc_POSTN",
    "Exc_SALL1_TDO2", "Exc_SATB2_TBR1", "Exc_SFTA3_NKX2-1",
    "Inh_EGFLAM", "Inh_EXPH5_TFAP2B", "Inh_OTX2_CASR", "Inh_PAX5_MCHR2",
    "Inh_PAX5_OTX1", "Inh_PAX5_POU3F1", "Inh_PENK", "Inh_PPP1R17_PCP2",
    "Inh_PPP1R1B_DRD3", "Inh_PVALB_SIX3", "Inh_SLC17A6_EBF3", "Inh_VIP_PROX1",
]

CELL_CLASS = {
    ct: ("Dopaminergic Neuron" if ct.startswith("Da_")
         else "Excitatory Neuron" if ct.startswith("Exc_")
         else "Inhibitory Neuron")
    for ct in CELL_TYPES
}

# ── Control-condition base proportions ────────────────────────────────────────
# Dopaminergic ~18%, Excitatory ~31%, Inhibitory ~51%
BASE_PROPS = {
    "Da_CALB1_CALCR":            0.020,
    "Da_CALB1_CRYM":             0.020,
    "Da_CALB1_TMEM200A_GAD2":    0.012,
    "Da_CALB1_TMEM200A_GEM":     0.015,
    "Da_CALB1_TMEM200A_TRHR":    0.010,
    "Da_SOX6_AGTR1":             0.020,
    "Da_SOX6_AGTR1_EYA4_NXPH1": 0.012,
    "Da_SOX6_EYA4":              0.030,   # most depleted in PD
    "Da_SOX6_GFRA2_PART1":       0.020,
    "Da_SOX6_NXPH1":             0.015,
    "Exc_5HT":                   0.028,
    "Exc_CYP2J2_CRH":            0.028,
    "Exc_HDC":                   0.024,
    "Exc_LHX9_FBN2":             0.034,
    "Exc_LMX1A_PITX2":           0.030,
    "Exc_NEUROD1":                0.028,
    "Exc_PAX5_F13A1":            0.027,
    "Exc_POSTN":                 0.025,
    "Exc_SALL1_TDO2":            0.027,
    "Exc_SATB2_TBR1":            0.025,
    "Exc_SFTA3_NKX2-1":          0.024,
    "Inh_EGFLAM":                0.042,
    "Inh_EXPH5_TFAP2B":          0.040,
    "Inh_OTX2_CASR":             0.037,
    "Inh_PAX5_MCHR2":            0.037,
    "Inh_PAX5_OTX1":             0.037,
    "Inh_PAX5_POU3F1":           0.034,
    "Inh_PENK":                  0.032,
    "Inh_PPP1R17_PCP2":          0.037,
    "Inh_PPP1R1B_DRD3":          0.037,
    "Inh_PVALB_SIX3":            0.032,
    "Inh_SLC17A6_EBF3":          0.032,
    "Inh_VIP_PROX1":             0.032,
}

BRAIN_BANKS = ["GTEX", "NethBB", "OHSU", "Sepulveda", "UKBB", "Umaryland"]
STUDIES     = ["Estiar", "Kamath"]


def _pd_props(base: dict, da_depletion: float = 0.45,
              sox6_eya4_extra: float = 0.30) -> np.ndarray:
    """Return PD-modified proportion vector (renormalised to 1)."""
    p = base.copy()
    for ct in p:
        if ct.startswith("Da_"):
            if ct == "Da_SOX6_EYA4":
                # Additional depletion on top of the general Da depletion
                p[ct] *= (1 - da_depletion) * (1 - sox6_eya4_extra)
            else:
                p[ct] *= (1 - da_depletion)
    total = sum(p.values())
    arr = np.array([p[ct] / total for ct in CELL_TYPES])
    return arr


def simulate_dataset(
    n_donors: int = 40,
    n_pd: int = 20,
    cells_per_sample: tuple = (150, 400),
    seed: int = 42,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    base_arr = np.array([BASE_PROPS[ct] for ct in CELL_TYPES])
    base_arr /= base_arr.sum()
    pd_arr = _pd_props(BASE_PROPS)

    rows = []

    for i in range(n_donors):
        donor_id    = f"D{i+1:03d}"
        is_pd       = i < n_pd
        case_ctrl   = "pd" if is_pd else "ctr"
        age         = float(rng.integers(58, 84) if is_pd else rng.integers(45, 80))
        sex         = rng.choice(["Male", "Female"], p=[0.60, 0.40])
        dapi_nurr   = rng.choice(["DAPI", "NURR"], p=[0.55, 0.45])
        brain_bank  = rng.choice(BRAIN_BANKS)
        study       = rng.choice(STUDIES)

        # Donors contribute 1-2 samples from SN, VTA, or RRF regions.
        # ~15 % of donors have only an RRF sample, breaking perfect collinearity
        # between region_SN and region_VTA across the cohort.
        region_draw = rng.random()
        if region_draw < 0.55:
            regions = [("SN", 1, 0), ("VTA", 0, 1)]   # both SN + VTA
        elif region_draw < 0.75:
            regions = [("SN", 1, 0)]                    # SN only
        elif region_draw < 0.90:
            regions = [("VTA", 0, 1)]                   # VTA only
        else:
            regions = [("RRF", 0, 0)]                   # RRF (neither)

        for region_name, region_SN, region_VTA in regions:
            sample_id = f"{donor_id}_{region_name}"
            n_cells   = int(rng.integers(*cells_per_sample))

            # Per-donor Dirichlet noise: concentration ∝ base proportion
            mean_props = pd_arr if is_pd else base_arr
            concentration = mean_props * 40          # lower = more variance
            donor_props = rng.dirichlet(concentration)
            cell_counts = rng.multinomial(n_cells, donor_props)

            for ct_idx, count in enumerate(cell_counts):
                ct = CELL_TYPES[ct_idx]
                for _ in range(count):
                    rows.append({
                        "sample_id":    sample_id,
                        "donor_id":     donor_id,
                        "case_control": case_ctrl,
                        "age":          age,
                        "sex":          sex,
                        "region_SN":    region_SN,
                        "region_VTA":   region_VTA,
                        "dapi_nurr":    dapi_nurr,
                        "brain_bank":   brain_bank,
                        "study":        study,
                        "cell_type":    ct,
                        "cell_class":   CELL_CLASS[ct],
                    })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Simulate snRNA-seq cell-composition data for crumblr testing."
    )
    parser.add_argument("--output",    required=True, help="Output CSV path")
    parser.add_argument("--n-donors",  type=int, default=40)
    parser.add_argument("--n-pd",      type=int, default=20)
    parser.add_argument("--min-cells", type=int, default=150)
    parser.add_argument("--max-cells", type=int, default=400)
    parser.add_argument("--seed",      type=int, default=42)
    args = parser.parse_args()

    df = simulate_dataset(
        n_donors=args.n_donors,
        n_pd=args.n_pd,
        cells_per_sample=(args.min_cells, args.max_cells),
        seed=args.seed,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    n_samples = df["sample_id"].nunique()
    n_donors  = df["donor_id"].nunique()
    size_mb   = out.stat().st_size / 1e6
    print(f"Simulated {len(df):,} cells | {n_donors} donors | {n_samples} samples")
    print(f"PD: {(df['case_control']=='pd').sum():,} cells | "
          f"CTR: {(df['case_control']=='ctr').sum():,} cells")
    print(f"Written to {out}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
