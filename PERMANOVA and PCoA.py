"""
PERMANOVA + PERMDISP + PCoA of defense module composition.
"""

import re
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2, f as f_dist

# =============================================================================
# User settings
# =============================================================================
FILE_PATH    = 'Matrix.csv'
PERMUTATIONS = 9999
SEED         = 42

# 'confidence' -> 95% confidence ellipse for the genus mean (matches manuscript)
# 'data'       -> 95% data / probability ellipse for the observations
ELLIPSE_TYPE = 'confidence'

# =============================================================================
# 1. Read CSV
# =============================================================================
df = pd.read_csv(FILE_PATH, sep=',', dtype=str).fillna('0')
df.rename(columns={df.columns[0]: 'Bac_name'}, inplace=True)

# Ensure unique column names
if df.columns.duplicated().any():
    print("Warning: duplicate column names detected - making them unique.")
    seen = {}
    new_cols = []
    for c in df.columns:
        if c in seen:
            seen[c] += 1
            new_cols.append(f"{c}.{seen[c]}")
        else:
            seen[c] = 0
            new_cols.append(c)
    df.columns = new_cols

module_cols = df.columns[1:].tolist()

for c in module_cols:
    df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)

print(f"Loaded {len(df)} genomes, {len(module_cols)} defense modules.")

# Sanity check for binary input
non_binary_counts = (~df[module_cols].isin([0, 1])).sum()
if (non_binary_counts > 0).any():
    print("\nWarning: values other than 0/1 detected. They will be treated as "
          "presence if > 0.")
    print(non_binary_counts[non_binary_counts > 0])

# =============================================================================
# 2. Genus assignment
# =============================================================================
GENUS_PREFIXES = [
    'Candidatus Epulonipiscioides',
    'Candidatus Epulonipiscium',
    'Candidatus Parepulonipiscium',
    'Candidatus Thiomargarita',
    'Epulopiscium',
    'Beggiatoa',
    'Thioploca',
    'Thiomargarita',
]

def assign_genus(name):
    name = str(name).strip()
    name_clean = re.sub(r'^uncultured\s+', '', name, flags=re.IGNORECASE)
    for prefix in GENUS_PREFIXES:
        if name_clean.startswith(prefix):
            return prefix
    m = re.search(r'(Candidatus\s+\w+)', name_clean)
    return m.group(1) if m else 'Other'

df['Genus'] = df['Bac_name'].apply(assign_genus)

genus_counts = df['Genus'].value_counts().sort_index()
print("\nGenus sample sizes:")
print(genus_counts)

if (genus_counts < 2).any():
    print("\nWarning: genera with fewer than 2 genomes present. "
          "PERMDISP may be unreliable for these groups:")
    print(genus_counts[genus_counts < 2])

# =============================================================================
# 3. Jaccard distance matrix
# =============================================================================
X_bin = (df[module_cols].values > 0).astype(np.uint8)
D_vec = pdist(X_bin, metric='jaccard')
D     = squareform(D_vec)

labels = df['Genus'].values
n      = len(labels)

# =============================================================================
# 4. PERMANOVA (Anderson 2001)
# =============================================================================
def permanova(D, labels, permutations=9999, seed=42):
    n = D.shape[0]
    labels = np.asarray(labels)
    D2 = D ** 2
    groups = np.unique(labels)
    a = len(groups)

    # Total sum of squares
    iu = np.triu_indices(n, k=1)
    SST = np.sum(D2[iu]) / n

    counts = Counter(labels)
    size_of_sample = np.array([counts[l] for l in labels], dtype=float)

    def ssw(lab, size_vec):
        same = (lab[:, None] == lab[None, :])
        return 0.5 * np.sum(D2 * same / size_vec[:, None])

    SSW = ssw(labels, size_of_sample)
    SSA = SST - SSW
    df_A = a - 1
    df_W = n - a

    F_obs = (SSA / df_A) / (SSW / df_W) if SSW > 0 else np.inf

    rng = np.random.default_rng(seed)
    F_null = np.empty(permutations)
    for k in range(permutations):
        perm = rng.permutation(n)
        lab_p = labels[perm]
        size_p = size_of_sample[perm]
        SSW_p = ssw(lab_p, size_p)
        SSA_p = SST - SSW_p
        F_null[k] = (SSA_p / df_A) / (SSW_p / df_W) if SSW_p > 0 else np.inf

    p_val = (np.sum(F_null >= F_obs) + 1) / (permutations + 1)

    # R^2
    R2 = SSA / SST if SST > 0 else np.nan
    return F_obs, p_val, R2, df_A, df_W, a

F_perm, p_perm, R2_perm, df_A, df_W, n_groups = permanova(
    D, labels, PERMUTATIONS, SEED
)

print("\n===== PERMANOVA =====")
print(f"Pseudo-F        = {F_perm:.4f}")
print(f"p-value         = {p_perm:.4f}   ({PERMUTATIONS} permutations)")
print(f"R^2             = {R2_perm:.4f}")
print(f"df (between)    = {df_A}")
print(f"df (within)     = {df_W}")
print(f"Number of groups = {n_groups}")

# =============================================================================
# 5. PCoA (classical MDS)
# =============================================================================
def pcoa(D):
    n = D.shape[0]
    A = -0.5 * D ** 2
    row_mean = A.mean(axis=1, keepdims=True)
    col_mean = A.mean(axis=0, keepdims=True)
    grand    = A.mean()
    G = A - row_mean - col_mean + grand

    eigvals, eigvecs = np.linalg.eigh(G)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pos_mask = eigvals > 1e-10
    eigvals_p = eigvals[pos_mask]
    eigvecs_p = eigvecs[:, pos_mask]

    coords = eigvecs_p * np.sqrt(eigvals_p)
    var_exp = eigvals_p / eigvals_p.sum()
    return coords, var_exp, eigvals

coords_full, var_exp_full, eigvals_all = pcoa(D)

n_pos = int(np.sum(eigvals_all >  1e-10))
n_neg = int(np.sum(eigvals_all < -1e-10))
print(f"\nPCoA eigenvalues: {n_pos} positive, {n_neg} negative")
if n_neg > 0:
    print("Note: negative eigenvalues present (Jaccard is not Euclidean). "
          "Variance percentages are based on positive eigenvalues only.")

coords = coords_full[:, :2]
var_exp = var_exp_full[:2] * 100.0

print(f"PC1 explains {var_exp[0]:.2f}% of variance, "
      f"PC2 explains {var_exp[1]:.2f}%")

# =============================================================================
# 6. PERMDISP (Anderson 2006)
# =============================================================================
def permdisp(coords, labels, permutations=9999, seed=42):
    n = coords.shape[0]
    labels = np.asarray(labels)
    groups = np.unique(labels)
    a = len(groups)

    # Distance to group centroid in full PCoA space
    dists = np.zeros(n)
    for g in groups:
        idx = np.where(labels == g)[0]
        if len(idx) == 0:
            continue
        centroid = coords[idx].mean(axis=0)
        dists[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)

    grand_mean = dists.mean()

    def anova_F(lab):
        SST = np.sum((dists - grand_mean) ** 2)
        SSW = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]
            if len(idx) == 0:
                continue
            gmean = dists[idx].mean()
            SSW += np.sum((dists[idx] - gmean) ** 2)
        SSA = SST - SSW
        df_A = a - 1
        df_W = n - a
        if SSW <= 0 or df_W <= 0:
            return np.inf
        return (SSA / df_A) / (SSW / df_W)

    F_obs = anova_F(labels)

    rng = np.random.default_rng(seed)
    F_null = np.empty(permutations)
    for k in range(permutations):
        F_null[k] = anova_F(rng.permutation(labels))

    p_val = (np.sum(F_null >= F_obs) + 1) / (permutations + 1)
    return F_obs, p_val, dists

F_disp, p_disp, dists_to_centroid = permdisp(
    coords_full, labels, PERMUTATIONS, SEED
)

print("\n===== PERMDISP =====")
print(f"F-stat          = {F_disp:.4f}")
print(f"p-value         = {p_disp:.4f}   ({PERMUTATIONS} permutations)")

print("\nInterpretation hints:")
if p_perm < 0.05 and p_disp >= 0.05:
    print(" - PERMANOVA significant + PERMDISP non-significant: "
          "strong support for differences in group centroids.")
elif p_perm < 0.05 and p_disp < 0.05:
    print(" - PERMANOVA significant + PERMDISP significant: "
          "dispersion differences may contribute; interpret PERMANOVA cautiously.")
elif p_perm >= 0.05:
    print(" - PERMANOVA non-significant: no evidence of centroid differences.")

# =============================================================================
# 7. Plot PCoA with ellipses
# =============================================================================
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 150})
sns.set_theme(style="whitegrid")

fig, ax = plt.subplots(figsize=(10, 8))

samples = pd.DataFrame(coords, columns=['PC1', 'PC2'], index=labels)
samples['Genus'] = labels

genus_list = sorted(np.unique(labels))
palette = dict(zip(genus_list, sns.color_palette("husl", len(genus_list))))

sns.scatterplot(
    data=samples, x='PC1', y='PC2',
    hue='Genus', palette=palette, hue_order=genus_list,
    s=60, alpha=0.85, edgecolor='black', linewidth=0.3, ax=ax
)

def add_ellipse(ax, data, color, ellipse_type='confidence'):
    if len(data) <= 2:
        return
    cov = np.cov(data[['PC1', 'PC2']].values.T)
    if np.any(~np.isfinite(cov)):
        return
    center = data[['PC1', 'PC2']].values.mean(axis=0)

    eigvals, eigvecs = np.linalg.eigh(cov)  # ascending
    major = eigvals[-1]
    minor = max(eigvals[0], 0.0)
    if major <= 0:
        return

    n_g = len(data)
    p   = 2

    if ellipse_type == 'data':
        scale = chi2.ppf(0.95, df=p)
    elif ellipse_type == 'confidence':
        if n_g <= p:
            return
        F_crit = f_dist.ppf(0.95, p, n_g - p)
        scale  = (p * (n_g - 1) / (n_g - p) * F_crit) / n_g
    else:
        raise ValueError("ellipse_type must be 'data' or 'confidence'")

    width  = 2.0 * np.sqrt(scale * major)
    height = 2.0 * np.sqrt(scale * minor)

    major_vec = eigvecs[:, -1]
    angle = np.degrees(np.arctan2(major_vec[1], major_vec[0]))

    ell = Ellipse(
        xy=center, width=width, height=height, angle=angle,
        edgecolor=color, facecolor='none',
        linewidth=1.5, linestyle='--'
    )
    ax.add_patch(ell)

for g in genus_list:
    sub = samples[samples['Genus'] == g]
    add_ellipse(ax, sub, palette[g], ellipse_type=ELLIPSE_TYPE)

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)', fontsize=12)
ax.set_title('PCoA of defense module composition (Jaccard distance)',
             fontsize=14, weight='bold')

ax.legend(
    title='Genus',
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    fontsize=9
)

plt.tight_layout()
plt.savefig('Defense_PCoA.svg', format='svg',
            bbox_inches='tight', facecolor='white')
plt.savefig('Defense_PCoA.png', format='png',
            dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

print("\nSaved: Defense_PCoA.svg and Defense_PCoA.png")
