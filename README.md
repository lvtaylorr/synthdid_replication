# Replication: Table 1 from Arkhangelsky et al. (2021)

This repository replicates Table 1 of **"Synthetic Difference-in-Differences"**
(Arkhangelsky, Athey, Hirshberg, Imbens, and Wager, *American Economic Review*, 2021).
The target result is the five treatment-effect estimates (SDID, SC, DID, DIFP, MC)
for the effect of California's Proposition 99 cigarette tax on per-capita cigarette
sales, along with placebo standard errors.

## Citation

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021).
Synthetic difference-in-differences. *American Economic Review*, 111(12), 4088–4118.
https://doi.org/10.1257/aer.20190159

## Data

The raw data (`input/california_prop99.csv`) are the California Proposition 99 panel
originally compiled by Abadie, Diamond, and Hainmueller (2010) and distributed with
the authors' replication package (OpenICPSR DOI: 10.3886/E146381V1).
The file is committed directly (< 25 KB). See `input/data_dictionary.md` for variable
descriptions.

## Prerequisites

- Python 3.9+ with `numpy`, `pandas`, `matplotlib`
- LaTeX (TeX Live or MacTeX) with `pdflatex` and `bibtex`
- GNU Make

Install Python dependencies:
```bash
pip install numpy pandas matplotlib
```

## Reproducing the paper

```bash
git clone https://github.com/YOUR-USERNAME/synthdid-sdid-paper.git
cd synthdid-sdid-paper
make
```

Or using the convenience wrapper:
```bash
bash run_all.sh
```

The paper compiles to `paper/paper.pdf`.

To rebuild from scratch:
```bash
make clean
make
```

## Pipeline

```
input/california_prop99.csv
        │
        ▼
code/preprocess.py  →  temp/panel.csv
        │
        ▼
code/analysis.py    →  output/tables/main_result.tex
                    →  output/figures/sdid_trends.pdf
        │
        ▼
paper/paper.tex     →  paper/paper.pdf
```

## Results

| Method | Replicated | Paper |
|--------|-----------|-------|
| SDID   | −15.6 (8.4)  | −15.6 (8.4)  |
| SC     | −19.5 (9.7)  | −19.6 (9.9)  |
| DID    | −27.3 (15.3) | −27.3 (17.7) |
| DIFP   | −11.1 (9.5)  | −11.1 (9.5)  |
| MC     | not replicated | −20.2 (11.5) |

SDID, SC, DID, and DIFP match the paper closely.
MC is not replicated because it requires the `MCPanel` R package (`mcnnm_cv`);
no equivalent Python implementation exists.

## Repository structure

```
synthdid-sdid-paper/
├── input/                  # Raw data (read-only)
│   ├── california_prop99.csv
│   └── data_dictionary.md
├── code/
│   ├── preprocess.py       # Load raw data → temp/panel.csv
│   └── analysis.py         # Estimate effects → output/
├── output/
│   ├── figures/
│   │   └── sdid_trends.pdf
│   └── tables/
│       └── main_result.tex
├── temp/                   # Intermediate files (gitignored)
├── paper/
│   ├── paper.tex
│   └── references.bib
├── Makefile
├── run_all.sh
├── proposal.md
└── README.md
```
