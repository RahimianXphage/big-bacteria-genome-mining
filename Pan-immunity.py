import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ----------------------------------------------------------------------
# 0. Settings
# ----------------------------------------------------------------------
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 9
rcParams['axes.linewidth'] = 0.8
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
rcParams['svg.fonttype'] = 'none'

base_dir = r"C:\Users\Qara\Desktop\Update 1\Ecothesis"

matrix_path = os.path.join(base_dir, "Matrix.csv")
bac_path    = os.path.join(base_dir, "bac data.xlsx")

np.random.seed(42)

N_PERM         = 999
N_ACCUM_ITER   = 500
MIN_PREVALENCE = 3

# ----------------------------------------------------------------------
# 1. Read and clean
# ----------------------------------------------------------------------
def clean_name(s: pd.Series) -> pd.Series:
    return (
        s.astype(str)
         .str.replace("\u00a0", " ", regex=False)
         .str.replace("\u200b", "",  regex=False)
         .str.replace("\u200c", "",  regex=False)
         .str.replace("\u200d", "",  regex=False)
         .str.strip()
    )

mat = pd.read_csv(matrix_path)
bac = pd.read_excel(bac_path)

mat["Bac name"]      = clean_name(mat["Bac name"])
bac["Bacteria name"] = clean_name(bac["Bacteria name"])

defense_cols_all = [c for c in mat.columns if c != "Bac name"]

M = mat.copy()
for c in defense_cols_all:
    M[c] = (M[c].astype(float) > 0).astype(int)

# ----------------------------------------------------------------------
# 2. Merge and verify metadata
# ----------------------------------------------------------------------
bac_slim = (
    bac[["Bacteria name", "Lifestyle"]]
    .rename(columns={"Bacteria name": "Bac name"})
)
df = M.merge(bac_slim, on="Bac name", how="left")

missing = int(df["Lifestyle"].isna().sum())
if missing > 0:
    print(f"[WARNING] {missing} genomes have no metadata match:")
    for n in df.loc[df["Lifestyle"].isna(), "Bac name"].tolist():
        print("   ", n)
    raise SystemExit(
        "Fix name mismatches before running the analysis. "
        "Do not silently drop or convert missing metadata."
    )
print(f"All {len(df)} genomes matched metadata.")

df["Lifestyle"] = df["Lifestyle"].astype(str).str.strip().str.lower()
df["Lifestyle"] = df["Lifestyle"].replace({
    "free-living": "Free-living", "free living": "Free-living",
    "freeliving":  "Free-living", "free":        "Free-living",
    "symbiotic":   "Symbiotic",   "symbiont":    "Symbiotic",
    "symbiosis":   "Symbiotic",   "sym":         "Symbiotic",
})
print("\nLifestyles:")
print(df["Lifestyle"].value_counts().to_string())

# ----------------------------------------------------------------------
# 3. Prevalence filter
# ----------------------------------------------------------------------
keep = [c for c in defense_cols_all if df[c].sum() >= MIN_PREVALENCE]
print(f"\nDefense modules retained (present in >= {MIN_PREVALENCE} genomes): "
      f"{len(keep)} of {len(defense_cols_all)}")

presence  = df[keep].to_numpy(dtype=np.int8)
lifestyle = df["Lifestyle"].to_numpy()

# ----------------------------------------------------------------------
# 4. Accumulation curves
# ----------------------------------------------------------------------
def accumulation_curve(sub_matrix, n_iter=N_ACCUM_ITER):
    n_genomes = sub_matrix.shape[0]
    if n_genomes < 2:
        return None, None
    curves = np.empty((n_iter, n_genomes), dtype=float)
    for it in range(n_iter):
        order   = np.random.permutation(n_genomes)
        running = np.maximum.accumulate(sub_matrix[order], axis=0)
        curves[it] = running.sum(axis=1)
    return curves.mean(axis=0), curves.std(axis=0)

print("\nComputing accumulation curves...")
mean_all,  sd_all  = accumulation_curve(presence)
mean_free, sd_free = accumulation_curve(presence[lifestyle == "Free-living"])
mean_symb, sd_symb = accumulation_curve(presence[lifestyle == "Symbiotic"])

n_all  = presence.shape[0]
n_free = int((lifestyle == "Free-living").sum())
n_symb = int((lifestyle == "Symbiotic").sum())

x_all  = np.arange(1, n_all  + 1)
x_free = np.arange(1, n_free + 1)
x_symb = np.arange(1, n_symb + 1)

max_n = max(n_all, n_free, n_symb)
def pad(arr, n):
    out = np.full(max_n, np.nan)
    out[:n] = arr
    return out

accum_df = pd.DataFrame({
    "n_genomes": x_all,
    "all_mean":  pad(mean_all,  n_all),
    "all_sd":    pad(sd_all,    n_all),
    "free_mean": pad(mean_free, n_free),
    "free_sd":   pad(sd_free,   n_free),
    "symb_mean": pad(mean_symb, n_symb),
    "symb_sd":   pad(sd_symb,   n_symb),
})
accum_df.to_csv(
    os.path.join(base_dir, "DefenseCooccurrence_accumulation.csv"),
    index=False,
)
print("Saved: DefenseCooccurrence_accumulation.csv")

# ----------------------------------------------------------------------
# 5. C-score with row/column-preserving null
# ----------------------------------------------------------------------
def c_score(matrix):
    """Average C-score (Stone & Roberts 1990).

    Equivalent closed form:
        D  = r[:, None] - M^T M,   r = column sums, D_ii = 0
        C  = trace(D @ D) / (n_cols * (n_cols - 1))
    """
    n_cols = matrix.shape[1]
    if n_cols < 2:
        return np.nan
    r = matrix.sum(axis=0, dtype=np.int64)
    S = matrix.T.astype(np.int64) @ matrix.astype(np.int64)
    D = r[:, None] - S
    return float(np.trace(D @ D)) / (n_cols * (n_cols - 1))


def swap_once(m, n_rows, n_cols):
    """One 2x2 checkerboard swap. Picks two rows first, then looks for a
    complementary 1/0 pair of columns. Mutates m in place. Returns True if
    a swap was made."""
    r1 = np.random.randint(n_rows)
    r2 = np.random.randint(n_rows)
    if r1 == r2:
        return False

    row1 = m[r1]
    row2 = m[r2]

    # Columns where row1 == 1 and row2 == 0
    cand_10 = np.flatnonzero((row1 == 1) & (row2 == 0))
    if cand_10.size == 0:
        return False
    # Columns where row1 == 0 and row2 == 1
    cand_01 = np.flatnonzero((row1 == 0) & (row2 == 1))
    if cand_01.size == 0:
        return False

    c1 = cand_10[np.random.randint(cand_10.size)]
    c2 = cand_01[np.random.randint(cand_01.size)]

    m[r1, c1] = 0; m[r1, c2] = 1
    m[r2, c1] = 1; m[r2, c2] = 0
    return True


def swap_randomize(matrix, n_requested):
    m = matrix.copy()
    n_rows, n_cols = m.shape
    n_cells = n_rows * n_cols
    max_attempts = max(n_requested * 20, n_cells * 100)

    accepted = 0
    attempts = 0
    while accepted < n_requested and attempts < max_attempts:
        if swap_once(m, n_rows, n_cols):
            accepted += 1
        attempts += 1
    return m, accepted, attempts


def c_score_test(matrix, n_perm=N_PERM):
    obs     = c_score(matrix)
    n_cells = matrix.size

    burn_in = 5 * n_cells
    between = 2 * n_cells

    m_chain, acc_burn, att_burn = swap_randomize(matrix, burn_in)

    null = np.empty(n_perm)
    total_accepted = acc_burn
    total_attempts = att_burn
    for p in range(n_perm):
        m_chain, acc, att = swap_randomize(m_chain, between)
        total_accepted += acc
        total_attempts += att
        null[p] = c_score(m_chain)

    mean_null = null.mean()
    sd_null   = null.std(ddof=1)
    ses       = (obs - mean_null) / sd_null if sd_null > 0 else np.nan
    p_one     = (1 + int(np.sum(null >= obs))) / (1 + n_perm)

    return {
        "observed_C":       obs,
        "null_mean":        mean_null,
        "null_sd":          sd_null,
        "SES":              ses,
        "p_one_sided":      p_one,
        "n_perm":           n_perm,
        "n_modules":        matrix.shape[1],
        "n_genomes":        matrix.shape[0],
        "swaps_requested":  burn_in + n_perm * between,
        "swaps_accepted":   total_accepted,
        "swaps_attempted":  total_attempts,
        "swap_accept_rate": (total_accepted / total_attempts
                             if total_attempts else np.nan),
    }


print(f"\nRunning C-score null model ({N_PERM} permutations)...")
print("  All genomes...")
res_all = c_score_test(presence)
print(f"    observed C = {res_all['observed_C']:.3f}")
print(f"    null mean  = {res_all['null_mean']:.3f} +/- {res_all['null_sd']:.3f}")
print(f"    SES        = {res_all['SES']:.2f}")
print(f"    p (1-sided)= {res_all['p_one_sided']:.4f}")
print(f"    swap accept rate = {res_all['swap_accept_rate']:.3f}")

results = {"All": res_all}

for ls in ["Free-living", "Symbiotic"]:
    sub = presence[lifestyle == ls]
    if sub.shape[0] >= 5 and sub.shape[1] >= 2:
        print(f"  {ls} (n = {sub.shape[0]})...")
        res = c_score_test(sub)
        print(f"    observed C = {res['observed_C']:.3f}")
        print(f"    null mean  = {res['null_mean']:.3f} +/- {res['null_sd']:.3f}")
        print(f"    SES        = {res['SES']:.2f}")
        print(f"    p (1-sided)= {res['p_one_sided']:.4f}")
        print(f"    swap accept rate = {res['swap_accept_rate']:.3f}")
        results[ls] = res

diag_rows = []
for k, v in results.items():
    diag_rows.append({
        "Group":            k,
        "observed_C":       v["observed_C"],
        "null_mean":        v["null_mean"],
        "null_sd":          v["null_sd"],
        "SES":              v["SES"],
        "p_one_sided":      v["p_one_sided"],
        "n_perm":           v["n_perm"],
        "n_genomes":        v["n_genomes"],
        "n_modules":        v["n_modules"],
        "swaps_requested":  v["swaps_requested"],
        "swaps_accepted":   v["swaps_accepted"],
        "swaps_attempted":  v["swaps_attempted"],
        "swap_accept_rate": v["swap_accept_rate"],
    })

cscore_df = pd.DataFrame(diag_rows).set_index("Group")
cscore_df.to_csv(os.path.join(base_dir, "DefenseCooccurrence_Cscore.csv"))
print("\nSaved: DefenseCooccurrence_Cscore.csv")

# ----------------------------------------------------------------------
# 6. Figure
# ----------------------------------------------------------------------
fig, axes = plt.subplots(
    1, 2, figsize=(9.0, 3.6),
    gridspec_kw={"width_ratios": [1.15, 1.0]},
)

ax = axes[0]
ax.plot(x_all,  mean_all,  color="#2c3e50", lw=1.6, label=f"All (n={n_all})")
ax.fill_between(x_all, mean_all - sd_all, mean_all + sd_all,
                color="#2c3e50", alpha=0.12)
ax.plot(x_free, mean_free, color="#3498db", lw=1.6,
        label=f"Free-living (n={n_free})")
ax.fill_between(x_free, mean_free - sd_free, mean_free + sd_free,
                color="#3498db", alpha=0.12)
ax.plot(x_symb, mean_symb, color="#e67e22", lw=1.6,
        label=f"Symbiotic (n={n_symb})")
ax.fill_between(x_symb, mean_symb - sd_symb, mean_symb + sd_symb,
                color="#e67e22", alpha=0.12)

ax.set_xlabel("Number of genomes sampled", fontsize=9)
ax.set_ylabel("Cumulative defense-module richness", fontsize=9)
ax.set_title("A. Defense-module accumulation (descriptive)",
             loc="left", fontsize=10, fontweight="bold", pad=8)
ax.legend(frameon=False, fontsize=8, loc="lower right")

ax = axes[1]
groups = list(results.keys())
obs    = [results[g]["observed_C"] for g in groups]
nm     = [results[g]["null_mean"]  for g in groups]
ns     = [results[g]["null_sd"]    for g in groups]
colors = ["#2c3e50", "#3498db", "#e67e22"][:len(groups)]
ypos   = np.arange(len(groups))[::-1]

for i, g in enumerate(groups):
    c = colors[i]
    ax.fill_betweenx([ypos[i] - 0.22, ypos[i] + 0.22],
                     nm[i] - 1.96 * ns[i], nm[i] + 1.96 * ns[i],
                     color=c, alpha=0.15)
    ax.plot(nm[i], ypos[i], marker="|", color=c, markersize=14, lw=0)
    ax.plot(obs[i], ypos[i], marker="o", color=c, markersize=8)
    ax.text(0.99, (ypos[i] + 0.55) / (len(groups) + 0.2),
            f"SES = {results[g]['SES']:.2f}, p = {results[g]['p_one_sided']:.3f}",
            transform=ax.transAxes, ha="right", fontsize=8, color="#333")

ax.set_yticks(ypos)
ax.set_yticklabels(groups)
ax.set_xlabel("C-score", fontsize=9)
ax.set_title("B. Checkerboard C-score vs row/column-preserving null",
             loc="left", fontsize=10, fontweight="bold", pad=8)

plt.tight_layout()
for ext in (".png", ".pdf", ".svg"):
    out = os.path.join(base_dir, f"DefenseCooccurrence_figure{ext}")
    fmt = ext.lstrip(".")
    plt.savefig(out, dpi=300 if fmt == "png" else None,
                bbox_inches="tight", format=fmt)
    print(f"Saved: {out}")

print("\nDone.")
