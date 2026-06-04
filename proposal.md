# Replication Proposal — GSE 552

## 1. Full citation

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W., & Wager, S. (2021).
Synthetic difference-in-differences. *American Economic Review*, 111(12), 4088–4118.
https://doi.org/10.1257/aer.20190159

## 2. Data source

California Proposition 99 panel dataset (Abadie, Diamond, and Hainmueller 2010),
distributed with the authors' replication package:
https://www.openicpsr.org/openicpsr/project/146381/version/V1/view
(DOI: 10.3886/E146381V1)

## 3. Result to replicate

I will replicate **Table 1** (all rows) from Arkhangelsky et al. (2021), which reports
treatment-effect estimates for the effect of California's Proposition 99 cigarette tax
on annual per-capita cigarette sales.
The numerical targets are:

| Method | Point estimate | Std. error |
|--------|---------------|-----------|
| SDID   | −15.6 packs   | 8.4       |
| SC     | −19.6 packs   | 9.9       |
| DID    | −27.3 packs   | 17.7      |
| DIFP   | −11.1 packs   | 9.5       |
| MC     | −20.2 packs   | 11.5      |

Note: MC requires the `MCPanel` R package and will be noted as a limitation
if a Python port is not feasible.

## 4. Planned toolchain

- **Language**: Python 3.11
- **Packages**: NumPy, Pandas, Matplotlib (standard; no non-standard packages)
- The Frank-Wolfe optimization algorithm will be implemented from scratch by
  translating the authors' R source (`R/solver.R`, `R/synthdid.R`), since no
  official Python package exists and the available third-party ports either
  fail to install or use different regularization parameters.
- LaTeX for the paper (`pdflatex` + `bibtex`).

## 5. Why this paper interests me

The SDID estimator elegantly bridges two dominant approaches in causal panel-data
analysis — synthetic control and difference-in-differences — while providing formal
asymptotic guarantees that neither method alone achieves.
I am interested in how the choice of estimator affects policy conclusions, and the
Proposition 99 application provides a concrete, well-known setting to evaluate that
question empirically.
