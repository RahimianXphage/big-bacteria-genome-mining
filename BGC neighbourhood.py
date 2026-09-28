from pathlib import Path
import re
import numpy as np
import pandas as pd

# ============================================================
# Configuration
# ============================================================
DATA = Path(r"D:\Projects\Gigantic\submit\Update 1\New island")

CONTIG_EDGE_MARGIN = 500       # bp: drop features within this of a contig end
CORE_TRIM_BP       = 3000      # bp: trim each side of antiSMASH region to approx core
N_PERM             = 10_000
RNG_SEED           = 42

SISTER_STRAINS = [
    ("Beggiatoa leptomitoformis strain D-401",
     "Beggiatoa leptomitoformis strain D-402"),
]

# ============================================================
# Helper functions
# ============================================================
def normalize_name(x):
    if pd.isna(x):
        return ""
    s = re.sub(r"\s+", " ", str(x).strip())
    return s.replace("'", "").replace('"', "")

def strip_contig_version(c):
    if pd.isna(c):
        return ""
    return re.sub(r"\.\d+$", "", str(c).strip())

def min_dist_to_core(def_start, def_stop, cores):
    """Minimum bp distance from a defense ORF to any BGC core on the same contig."""
    if len(cores) == 0:
        return np.nan
    dists = []
    for cs, ce in cores:
        if def_stop < cs:
            dists.append(cs - def_stop)
        elif def_start > ce:
            dists.append(def_start - ce)
        else:
            dists.append(0)
    return min(dists)

# ============================================================
# 1. Load inputs
# ============================================================
def_ = pd.read_excel(DATA / "Defense.xlsx")
bgc  = pd.read_excel(DATA / "BGC.xlsx")
def_.columns = [c.strip() for c in def_.columns]
bgc.columns  = [c.strip() for c in bgc.columns]

def_ = def_.rename(columns={
    "bacteria_name":       "genome",
    "contig_name":         "contig",
    "defense_system_name": "defense",
})
bgc = bgc.rename(columns={
    "Bacteria name": "genome",
    "Contig":        "contig",
    "Region":        "region",
    "Type":          "bcg_type",
    "From":          "region_start",
    "To":            "region_end",
})

bgc  = bgc[bgc["region"].astype(str).str.contains("Region", na=False)].copy()
bgc  = bgc.dropna(subset=["region_start", "region_end"]).copy()
def_ = def_.dropna(subset=["start", "stop"]).copy()

bgc["region_start"] = bgc["region_start"].astype(int)
bgc["region_end"]   = bgc["region_end"].astype(int)
def_["start"]       = def_["start"].astype(int)
def_["stop"]        = def_["stop"].astype(int)

def_["genome_norm"] = def_["genome"].apply(normalize_name)
def_["contig_norm"] = def_["contig"].apply(strip_contig_version)
def_["key"]         = list(zip(def_["genome_norm"], def_["contig_norm"]))

bgc["genome_norm"]  = bgc["genome"].apply(normalize_name)
bgc["contig_norm"]  = bgc["contig"].apply(strip_contig_version)
bgc["key"]          = list(zip(bgc["genome_norm"], bgc["contig_norm"]))

print(f"Raw defense rows : {len(def_)}")
print(f"Raw BGC regions  : {len(bgc)}")

# ============================================================
# 2. Collapse sister strains
# ============================================================
for keep, drop in SISTER_STRAINS:
    k = normalize_name(keep)
    d = normalize_name(drop)
    n_before = len(def_) + len(bgc)
    def_ = def_[def_["genome_norm"] != d].copy()
    bgc  = bgc[bgc["genome_norm"]  != d].copy()
    print(f"Collapsed {drop} -> {keep}: removed "
          f"{n_before - (len(def_) + len(bgc))} rows")

# ============================================================
# 3. Contig-edge filter
# ============================================================
def_before = len(def_)
bgc_before = len(bgc)
def_ = def_[def_["start"] >= CONTIG_EDGE_MARGIN].copy()
bgc  = bgc[bgc["region_start"] >= CONTIG_EDGE_MARGIN].copy()

print(f"\nContig-edge filter (>= {CONTIG_EDGE_MARGIN} bp):")
print(f"  Defense kept {len(def_)} / {def_before}")
print(f"  BGC     kept {len(bgc)} / {bgc_before}")

# ============================================================
# 4. Approximate BGC core biosynthetic span
# ============================================================
bgc["core_start"] = bgc["region_start"] + CORE_TRIM_BP
bgc["core_end"]   = bgc["region_end"]   - CORE_TRIM_BP
bgc = bgc[bgc["core_end"] > bgc["core_start"]].copy()

print(f"\nCore-span approximation (trim {CORE_TRIM_BP/1000:.0f} kb per side):")
print(f"  Regions with valid core: {len(bgc)}")

# ============================================================
# 5. Classify defense ORFs
# ============================================================
bgc_by_key = {k: v for k, v in bgc.groupby("key")}

def classify(row):
    regs = bgc_by_key.get(row["key"])
    if regs is None or len(regs) == 0:
        return "none", None, None
    hits = []
    for _, r in regs.iterrows():
        if row["start"] >= r["core_start"] and row["stop"] <= r["core_end"]:
            hits.append(("integrated", r["region"], r["bcg_type"]))
        elif row["start"] >= r["region_start"] and row["stop"] <= r["region_end"]:
            hits.append(("adjacent", r["region"], r["bcg_type"]))
    if not hits:
        return "none", None, None
    integ = [h for h in hits if h[0] == "integrated"]
    return integ[0] if integ else hits[0]

rows = []
for _, d in def_.iterrows():
    status, region, bcg_type = classify(d)
    rows.append({
        "Genome":     d["genome"],
        "Contig":     d["contig"],
        "Defense":    d["defense"],
        "Defense_start": d["start"],
        "Defense_stop":  d["stop"],
        "Status":     status,
        "BGC_region": region,
        "BGC_type":   bcg_type,
    })
reannotated = pd.DataFrame(rows)

n_int = int((reannotated["Status"] == "integrated").sum())
n_adj = int((reannotated["Status"] == "adjacent").sum())
n_non = int((reannotated["Status"] == "none").sum())

print(f"\nDefense ORF classification:")
print(f"  Integrated : {n_int}")
print(f"  Adjacent   : {n_adj}")
print(f"  None       : {n_non}")

# ============================================================
# 6. Keep only the island records
# ============================================================
islands = reannotated[reannotated["Status"].isin(["integrated", "adjacent"])].copy()
islands = islands.sort_values(
    ["Genome", "Contig", "Defense_start", "Status"]
).reset_index(drop=True)

print(f"\nIsland records written to output: {len(islands)}")

# ============================================================
# 7. Null model: within-contig permutation
# ============================================================
# Contig max coordinate = approximate contig length
contig_max = {}
for _, r in def_.iterrows():
    k = r["key"]
    contig_max[k] = max(contig_max.get(k, 0), int(r["stop"]))
for _, r in bgc.iterrows():
    k = r["key"]
    contig_max[k] = max(contig_max.get(k, 0), int(r["region_end"]))

shared_keys = sorted(set(def_["key"]) & set(bgc["key"]))
def_count   = (def_[def_["key"].isin(shared_keys)]
               .groupby("key").size().to_dict())
cores_by_key = {k: bgc[bgc["key"] == k][["core_start", "core_end"]].values
                for k in shared_keys}

print(f"\nShared contigs for null model: {len(shared_keys)}")
print(f"Total defenses on those contigs: {sum(def_count.values())}")

rng = np.random.default_rng(RNG_SEED)
null_counts = np.zeros(N_PERM, dtype=int)

for b in range(N_PERM):
    total = 0
    for k in shared_keys:
        n_def = def_count.get(k, 0)
        if n_def == 0:
            continue
        L = contig_max.get(k, 0)
        if L < 2 * CONTIG_EDGE_MARGIN + 1000:
            continue
        cores = cores_by_key.get(k)
        if cores is None or len(cores) == 0:
            continue
        pos = rng.integers(CONTIG_EDGE_MARGIN,
                           L - CONTIG_EDGE_MARGIN,
                           size=n_def)
        inside = np.zeros(n_def, dtype=bool)
        for cs, ce in cores:
            if ce <= cs:
                continue
            inside |= (pos >= cs) & (pos <= ce)
        total += int(inside.sum())
    null_counts[b] = total

obs  = n_int
pval = (1 + int((null_counts >= obs).sum())) / (N_PERM + 1)
oer  = obs / null_counts.mean() if null_counts.mean() else np.nan

print(f"\nNesting permutation test:")
print(f"  Observed integrated : {obs}")
print(f"  Null mean           : {null_counts.mean():.2f}")
print(f"  Null SD             : {null_counts.std():.2f}")
print(f"  Null 95th percentile: {np.percentile(null_counts, 95):.1f}")
print(f"  Observed / Expected : {oer:.3f}")
print(f"  Empirical p-value   : {pval:.4f}")

# ============================================================
# 8. Distance-sensitivity analysis (for reviewer's "50 kb" critique)
# ============================================================
def_lengths  = {k: (grp["stop"] - grp["start"]).values
                for k, grp in def_.groupby("key")}

# Observed distances
obs_dists = []
for _, r in def_.iterrows():
    cores = cores_by_key.get(r["key"])
    if cores is None or len(cores) == 0:
        continue
    d = min_dist_to_core(r["start"], r["stop"], cores)
    if not np.isnan(d):
        obs_dists.append(d)
obs_dists = np.array(obs_dists)

THRESHOLDS = [10_000, 25_000, 50_000, 100_000]
null_thresh = {X: np.zeros(N_PERM, dtype=int) for X in THRESHOLDS}

for b in range(N_PERM):
    for k in shared_keys:
        n_def = def_count.get(k, 0)
        if n_def == 0:
            continue
        L = contig_max.get(k, 0)
        if L < 2 * CONTIG_EDGE_MARGIN + 1000:
            continue
        cores = cores_by_key.get(k)
        if cores is None or len(cores) == 0:
            continue
        lens = def_lengths.get(k, np.array([500]))
        for _ in range(n_def):
            length = int(rng.choice(lens))
            start  = int(rng.integers(CONTIG_EDGE_MARGIN,
                                      L - CONTIG_EDGE_MARGIN - length))
            stop   = start + length
            d = min_dist_to_core(start, stop, cores)
            if np.isnan(d):
                continue
            for X in THRESHOLDS:
                if d <= X:
                    null_thresh[X][b] += 1

dist_rows = []
for X in THRESHOLDS:
    obs_n    = int((obs_dists <= X).sum())
    exp_n    = float(null_thresh[X].mean())
    oer_t    = obs_n / exp_n if exp_n else np.nan
    p_t      = (1 + int((null_thresh[X] >= obs_n).sum())) / (N_PERM + 1)
    dist_rows.append({
        "Threshold_kb": X // 1000,
        "Observed":     obs_n,
        "Expected":     round(exp_n, 2),
        "O/E":          round(oer_t, 3),
        "p":            round(p_t, 4),
    })
    print(f"  ≤{X//1000:>3d} kb  observed = {obs_n:>4d}  "
          f"expected = {exp_n:7.2f}  O/E = {oer_t:5.2f}  p = {p_t:.4f}")

dist_sens = pd.DataFrame(dist_rows)

# ============================================================
# 9. Build summary sheet
# ============================================================
n_islands_unique = islands[["Genome", "BGC_region"]].drop_duplicates().shape[0]

summary_df = pd.DataFrame({
    "Parameter": [
        "Contig-edge margin (bp)",
        "Core trim (kb)",
        "Permutations",
        "Genomes with integrated defense ORFs",
        "Defense ORFs classified as integrated",
        "Defense ORFs classified as adjacent",
        "Total island records in output",
        "Shared contigs in null model",
        "Nesting null mean",
        "Nesting null SD",
        "Nesting null 95th percentile",
        "Nesting observed / expected",
        "Nesting empirical p-value",
    ],
    "Value": [
        CONTIG_EDGE_MARGIN,
        CORE_TRIM_BP // 1000,
        N_PERM,
        islands[islands["Status"] == "integrated"]["Genome"].nunique(),
        n_int,
        n_adj,
        len(islands),
        len(shared_keys),
        round(float(null_counts.mean()), 2),
        round(float(null_counts.std()), 2),
        round(float(np.percentile(null_counts, 95)), 1),
        round(float(oer), 3),
        round(float(pval), 4),
    ],
})

# ============================================================
# 10. Save everything to one workbook
# ============================================================
out_xlsx = DATA / "islands_results.xlsx"

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
    islands.to_excel(w, sheet_name="Islands", index=False)
    bgc[["genome", "contig", "region", "bcg_type",
         "region_start", "region_end", "core_start", "core_end"]].to_excel(
        w, sheet_name="BGC_cores", index=False)
    pd.DataFrame({"null_integrated": null_counts}).to_excel(
        w, sheet_name="Null_distribution", index=False)
    dist_sens.to_excel(w, sheet_name="Distance_sensitivity", index=False)
    summary_df.to_excel(w, sheet_name="Summary", index=False)

print(f"\nSaved: {out_xlsx}")
print(f"Sheets: Islands | BGC_cores | Null_distribution | "
      f"Distance_sensitivity | Summary")
