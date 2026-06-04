"""
code/analysis.py

Reads temp/panel.csv (written by code/preprocess.py), estimates SDID, SC,
DID, and DIFP treatment effects, and writes:
  - output/tables/main_result.tex   (LaTeX table fragment)
  - output/figures/sdid_trends.pdf  (parallel-trends plot)

Run from the repository root:
    python code/analysis.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# ── Paths (all relative to repo root) ────────────────────────────────────────
PANEL_FILE  = "temp/panel.csv"
TABLE_OUT   = "output/tables/main_result.tex"
FIGURE_OUT  = "output/figures/sdid_trends.pdf"

os.makedirs("output/tables",  exist_ok=True)
os.makedirs("output/figures", exist_ok=True)

# ── Load panel ────────────────────────────────────────────────────────────────
panel = pd.read_csv(PANEL_FILE, index_col=0)
panel.columns = panel.columns.astype(int)

states = list(panel.index)
years  = list(panel.columns)

# California is last (ordered by preprocess.py)
N0 = len(states) - 1       # 38 control states
T0 = len([y for y in years if y < 1989])  # 19 pre-treatment years
Y  = panel.values.astype(float)           # shape (39, 31)
N, T = Y.shape
N1 = N - N0    # 1
T1 = T - T0    # 12


# ── Core helpers ──────────────────────────────────────────────────────────────

def collapsed_form(Y, N0, T0):
    top_left  = Y[:N0, :T0]
    top_right = Y[:N0, T0:].mean(axis=1, keepdims=True)
    bot_left  = Y[N0:, :T0].mean(axis=0, keepdims=True)
    bot_right = np.array([[Y[N0:, T0:].mean()]])
    return np.block([[top_left, top_right], [bot_left, bot_right]])


def fw_step(A, x, b, eta):
    Ax        = A @ x
    half_grad = (Ax - b) @ A + eta * x
    i         = np.argmin(half_grad)
    d_x       = -x.copy()
    d_x[i]   += 1.0
    if np.all(d_x == 0):
        return x
    d_err = A[:, i] - Ax
    step  = -(half_grad @ d_x) / (np.sum(d_err ** 2) + eta * np.sum(d_x ** 2))
    step  = float(np.clip(step, 0, 1))
    return x + step * d_x


def sc_weight_fw(Y_in, zeta, intercept=True, lam_init=None,
                 min_decrease=1e-3, max_iter=1000):
    N_obs, T_plus1 = Y_in.shape
    T_pred = T_plus1 - 1
    lam    = np.ones(T_pred) / T_pred if lam_init is None else lam_init.copy()
    Yw     = Y_in - Y_in.mean(axis=0) if intercept else Y_in.copy()
    A, b   = Yw[:, :T_pred], Yw[:, T_pred]
    eta    = N_obs * zeta ** 2
    prev_val = np.inf
    for t in range(max_iter):
        lam = fw_step(A, lam, b, eta)
        err = Yw @ np.append(lam, -1.0)
        val = zeta ** 2 * np.sum(lam ** 2) + np.sum(err ** 2) / N_obs
        if t >= 1 and prev_val - val < min_decrease ** 2:
            break
        prev_val = val
    return lam


def sparsify(v):
    v = v.copy()
    v[v <= v.max() / 4] = 0.0
    return v / v.sum()


def estimate_weights(Y, N0, T0, zeta_omega, zeta_lambda, noise_level,
                     omega_intercept=True, lambda_intercept=True,
                     max_iter=10_000, max_iter_pre_sparsify=100):
    min_dec  = 1e-5 * noise_level
    Yc       = collapsed_form(Y, N0, T0)
    lam_input = Yc[:N0, :]
    lam0 = sc_weight_fw(lam_input, zeta_lambda, intercept=lambda_intercept,
                        min_decrease=min_dec, max_iter=max_iter_pre_sparsify)
    lam  = sc_weight_fw(lam_input, zeta_lambda, intercept=lambda_intercept,
                        lam_init=sparsify(lam0), min_decrease=min_dec,
                        max_iter=max_iter)
    om_input = Yc[:, :T0].T
    om0 = sc_weight_fw(om_input, zeta_omega, intercept=omega_intercept,
                       min_decrease=min_dec, max_iter=max_iter_pre_sparsify)
    om  = sc_weight_fw(om_input, zeta_omega, intercept=omega_intercept,
                       lam_init=sparsify(om0), min_decrease=min_dec,
                       max_iter=max_iter)
    return lam, om


def point_estimate(Y, N0, T0, omega, lam):
    N1_   = Y.shape[0] - N0
    T1_   = Y.shape[1] - T0
    left  = np.append(-omega, np.ones(N1_) / N1_)
    right = np.append(-lam,   np.ones(T1_) / T1_)
    return float(left @ Y @ right)


# ── Noise level and regularization ───────────────────────────────────────────
noise_level = float(np.std(np.diff(Y[:N0, :T0], axis=1), ddof=1))
eta_omega   = (N1 * T1) ** 0.25
eta_lambda  = 1e-6
zeta_omega  = eta_omega  * noise_level
zeta_lambda = eta_lambda * noise_level
zeta_sc     = 1e-6 * noise_level


# ── SDID ──────────────────────────────────────────────────────────────────────
lam_sdid, om_sdid = estimate_weights(Y, N0, T0, zeta_omega, zeta_lambda,
                                      noise_level)
tau_sdid = point_estimate(Y, N0, T0, om_sdid, lam_sdid)

# ── DID ───────────────────────────────────────────────────────────────────────
tau_did = point_estimate(Y, N0, T0,
                         omega=np.ones(N0) / N0,
                         lam=np.ones(T0) / T0)

# ── SC ────────────────────────────────────────────────────────────────────────
Yc_sc       = collapsed_form(Y, N0, T0)
om_sc_pre   = sc_weight_fw(Yc_sc[:, :T0].T, zeta_sc, intercept=False,
                            max_iter=100)
om_sc       = sc_weight_fw(Yc_sc[:, :T0].T, zeta_sc, intercept=False,
                            lam_init=sparsify(om_sc_pre))
tau_sc      = point_estimate(Y, N0, T0, om_sc, np.zeros(T0))

# ── DIFP (de-meaned SC, Doudchenko & Imbens 2016) ─────────────────────────────
Yc_difp       = collapsed_form(Y, N0, T0)
om_difp_pre   = sc_weight_fw(Yc_difp[:, :T0].T, zeta_sc, intercept=True,
                              max_iter=100)
om_difp       = sc_weight_fw(Yc_difp[:, :T0].T, zeta_sc, intercept=True,
                              lam_init=sparsify(om_difp_pre))
tau_difp      = point_estimate(Y, N0, T0, om_difp, np.ones(T0) / T0)


# ── Placebo standard errors (Algorithm 4) ────────────────────────────────────
def placebo_se(Y, N0, T0, estimator="sdid", replications=200, seed=42):
    rng = np.random.default_rng(seed)
    N1_ = Y.shape[0] - N0
    placebo_ests = []
    for _ in range(replications):
        perm       = rng.permutation(N0)
        ctrl_idx   = perm[:N0 - N1_]
        pseudo_idx = perm[N0 - N1_:]
        Y_b  = np.vstack([Y[ctrl_idx], Y[pseudo_idx]])
        N0_b = len(ctrl_idx)
        if estimator == "sdid":
            lam_b, om_b = estimate_weights(Y_b, N0_b, T0, zeta_omega,
                                           zeta_lambda, noise_level,
                                           max_iter=1000)
            est = point_estimate(Y_b, N0_b, T0, om_b, lam_b)
        elif estimator == "did":
            est = point_estimate(Y_b, N0_b, T0,
                                 np.ones(N0_b) / N0_b,
                                 np.ones(T0) / T0)
        elif estimator == "sc":
            Yc_b = collapsed_form(Y_b, N0_b, T0)
            om_b0 = sc_weight_fw(Yc_b[:, :T0].T, zeta_sc,
                                  intercept=False, max_iter=100)
            om_b  = sc_weight_fw(Yc_b[:, :T0].T, zeta_sc,
                                  intercept=False, lam_init=sparsify(om_b0))
            est   = point_estimate(Y_b, N0_b, T0, om_b, np.zeros(T0))
        elif estimator == "difp":
            Yc_b = collapsed_form(Y_b, N0_b, T0)
            om_b0 = sc_weight_fw(Yc_b[:, :T0].T, zeta_sc,
                                  intercept=True, max_iter=100)
            om_b  = sc_weight_fw(Yc_b[:, :T0].T, zeta_sc,
                                  intercept=True, lam_init=sparsify(om_b0))
            est   = point_estimate(Y_b, N0_b, T0, om_b, np.ones(T0) / T0)
        placebo_ests.append(est)
    ests = np.array(placebo_ests)
    R    = len(ests)
    return float(np.sqrt((R - 1) / R) * np.std(ests, ddof=1))


se_sdid = placebo_se(Y, N0, T0, estimator="sdid", replications=200)
se_did  = placebo_se(Y, N0, T0, estimator="did",  replications=200)
se_sc   = placebo_se(Y, N0, T0, estimator="sc",   replications=200)
se_difp = placebo_se(Y, N0, T0, estimator="difp", replications=200)


# ── Write LaTeX table ─────────────────────────────────────────────────────────
table = r"""\begin{table}[h]
\centering
\caption{Effect of Proposition 99 on California Per-Capita Cigarette Sales}
\label{tab:main}
\small
\begin{tabular}{lccc}
\toprule
Method & Replication & Paper & Std.\ Error \\
\midrule
""" + \
f"Synthetic DiD (SDID)            & ${tau_sdid:.1f}$ & $-15.6$ & $({se_sdid:.1f})$ \\\\\n" + \
f"Synthetic Control (SC)          & ${tau_sc:.1f}$ & $-19.6$ & $({se_sc:.1f})$ \\\\\n" + \
f"Difference-in-Differences (DID) & ${tau_did:.1f}$ & $-27.3$ & $({se_did:.1f})$ \\\\\n" + \
f"De-meaned SC (DIFP)             & ${tau_difp:.1f}$ & $-11.1$ & $({se_difp:.1f})$ \\\\\n" + \
r"""Matrix Completion (MC)          & ---     & $-20.2$ & $(11.5)$ \\
\midrule
\multicolumn{4}{l}{\footnotesize Notes: Outcome is per-capita cigarette sales in packs per year.} \\
\multicolumn{4}{l}{\footnotesize Sample: 39 U.S.\ states, 1970--2000. Treatment: California post-1989.} \\
\multicolumn{4}{l}{\footnotesize Standard errors via placebo method (200 replications).} \\
\multicolumn{4}{l}{\footnotesize MC not replicated: depends on MCPanel R package (see text).} \\
\bottomrule
\end{tabular}
\end{table}
"""

with open(TABLE_OUT, "w") as f:
    f.write(table)


# ── Write parallel-trends figure ──────────────────────────────────────────────
# Synthetic California = SDID-weighted average of 38 donor states
synth_ca = Y[:N0, :].T @ om_sdid   # shape (31,)
actual_ca = Y[-1, :]                 # shape (31,)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(years, actual_ca, color="black",    linewidth=2,   label="California (actual)")
ax.plot(years, synth_ca,  color="steelblue", linewidth=2,
        linestyle="--", label="Synthetic California (SDID weights)")
ax.axvline(x=1989, color="gray", linestyle=":", linewidth=1.5, label="Proposition 99 (1989)")
ax.set_xlabel("Year")
ax.set_ylabel("Per-capita cigarette sales (packs/year)")
ax.set_title("California vs. Synthetic California, 1970–2000")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURE_OUT, dpi=300)
plt.close(fig)


# ── Console summary ───────────────────────────────────────────────────────────
print("=" * 55)
print("Table 1 Replication: Effect of Proposition 99")
print("Outcome: Per-capita cigarette sales (packs/year)")
print("=" * 55)
print(f"{'Method':<8}  {'Replicated':>12}  {'Paper':>8}  {'SE':>8}")
print("-" * 45)
print(f"{'SDID':<8}  {tau_sdid:>12.1f}  {-15.6:>8.1f}  ({se_sdid:.1f})")
print(f"{'SC':<8}  {tau_sc:>12.1f}  {-19.6:>8.1f}  ({se_sc:.1f})")
print(f"{'DID':<8}  {tau_did:>12.1f}  {-27.3:>8.1f}  ({se_did:.1f})")
print(f"{'DIFP':<8}  {tau_difp:>12.1f}  {-11.1:>8.1f}  ({se_difp:.1f})")
print(f"{'MC':<8}  {'---':>12}  {-20.2:>8.1f}  (not replicated)")
print("=" * 55)
print(f"Table written to : {TABLE_OUT}")
print(f"Figure written to: {FIGURE_OUT}")
