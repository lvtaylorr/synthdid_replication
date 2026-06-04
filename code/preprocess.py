"""
code/preprocess.py

Reads raw data from input/, cleans it, and writes an analysis-ready
panel matrix to temp/panel.csv.

Run from the repository root:
    python code/preprocess.py
"""

import numpy as np
import pandas as pd
import os

# ── Paths (all relative to repo root) ────────────────────────────────────────
INPUT_FILE  = "input/california_prop99.csv"
OUTPUT_FILE = "temp/panel.csv"

os.makedirs("temp", exist_ok=True)

# ── Load raw data ─────────────────────────────────────────────────────────────
# Note: the CSV uses a semicolon delimiter, not a comma.
df = pd.read_csv(INPUT_FILE, sep=";")

# Basic integrity checks
required_cols = {"State", "Year", "PacksPerCapita", "treated"}
assert required_cols.issubset(df.columns), f"Missing columns: {required_cols - set(df.columns)}"
assert df["PacksPerCapita"].notna().all(), "Missing values in PacksPerCapita"

# ── Build ordered panel: control states first, treated state (California) last ─
treat_by_state   = df.groupby("State")["treated"].max()
control_states   = sorted(treat_by_state[treat_by_state == 0].index.tolist())
treated_states   = treat_by_state[treat_by_state == 1].index.tolist()
ordered_states   = control_states + treated_states

# Pivot to wide matrix: rows = states (ordered), columns = years (sorted)
panel = df.pivot(index="State", columns="Year", values="PacksPerCapita")
panel = panel.loc[ordered_states]            # enforce state order
panel = panel[sorted(panel.columns)]         # enforce year order (1970–2000)

# Write to temp/ as CSV (index = State names, columns = years)
panel.to_csv(OUTPUT_FILE)

# ── Console summary ───────────────────────────────────────────────────────────
N, T    = panel.shape
N0      = len(control_states)
T0      = len([c for c in panel.columns if c < 1989])
years   = list(panel.columns)
outcome = df["PacksPerCapita"]

print("=" * 50)
print("Preprocessing complete")
print("=" * 50)
print(f"  Input file     : {INPUT_FILE}")
print(f"  Output file    : {OUTPUT_FILE}")
print(f"  States (N)     : {N}  ({N0} control + {N - N0} treated)")
print(f"  Years  (T)     : {T}  ({years[0]}–{years[-1]})")
print(f"  Pre-treatment  : T0 = {T0} ({years[0]}–{years[T0-1]})")
print(f"  Post-treatment : T1 = {T - T0} ({years[T0]}–{years[-1]})")
print(f"  Observations   : {N * T}")
print(f"  PacksPerCapita : mean = {outcome.mean():.1f}, "
      f"std = {outcome.std():.1f}, "
      f"min = {outcome.min():.1f}, "
      f"max = {outcome.max():.1f}")
print("=" * 50)
