from __future__ import annotations
import math
from typing import List, Dict

def find_spcas9_sites(seq: str, pam: str = "NGG") -> List[tuple[int,str,str]]:
    """
    Simple SpCas9 finder: 20nt guide immediately 5' of PAM (NGG).
    Returns list of tuples (position, guide, pam_trimer).
    """
    seq = seq.upper()
    hits = []
    for i in range(len(seq) - 23):
        guide = seq[i:i+20]
        pam_tri = seq[i+20:i+23]
        if len(guide) == 20 and len(pam_tri) == 3 and pam_tri.endswith("GG"):
            hits.append((i, guide, pam_tri))
    return hits

def gc_percent(s: str) -> float:
    if not s: return 0.0
    g = s.count("G") + s.count("C")
    return 100.0 * g / len(s)

def cpg_density(s: str) -> float:
    # count “CG” per 20bp (normalized so typical values are ~0–0.3)
    if not s: return 0.0
    n = sum(1 for i in range(len(s)-1) if s[i:i+2] == "CG")
    return round(n / len(s), 3)

def codon_bias_index(s: str) -> float:
    """
    Lightweight proxy: ratio of GC di-nucleotides to AT di-nucleotides in the guide.
    Returns ~0.8–1.2 typically; closer to 1 is neutral.
    """
    if len(s) < 2: return 1.0
    gc_pairs = sum(1 for i in range(len(s)-1) if s[i] in "GC" and s[i+1] in "GC")
    at_pairs = sum(1 for i in range(len(s)-1) if s[i] in "AT" and s[i+1] in "AT")
    return round((gc_pairs + 1) / (at_pairs + 1), 2)

def heuristic_offtarget_score(guide: str, cpg: float) -> float:
    """
    0..100; higher is better (lower off-target likelihood).
    Simple: penalize extreme GC and high CpG.
    """
    gc = gc_percent(guide) / 100.0
    gc_pen = abs(0.5 - gc) * 100  # deviation from 50% GC
    cpg_pen = cpg * 100 * 0.7     # weight CpG
    score = max(0.0, 100.0 - (gc_pen + cpg_pen))
    return round(score, 2)

def score_guides(seq: str, max_rows: int = 10) -> List[Dict]:
    guides = find_spcas9_sites(seq)
    rows = []
    for (pos, g, pam) in guides[:200]:  # cap for speed
        gc = round(gc_percent(g), 2)
        cpg = cpg_density(g)
        cb = codon_bias_index(g)
        ot = heuristic_offtarget_score(g, cpg)
        rows.append({
            "Position": pos,
            "Guide RNA": g,
            "PAM": pam,
            "GC Score": gc,
            "CpG": cpg,
            "Codon Bias": cb,
            "Off-target (↑=better)": ot
        })
    # sort by off-target desc, then GC centered near 50.
    rows.sort(key=lambda r: (r["Off-target (↑=better)"], -abs(50-r["GC Score"])), reverse=True)
    return rows[:max_rows]

def simulated_cut_preview(seq: str, first_hit_pos: int | None) -> str:
    if first_hit_pos is None or first_hit_pos < 10:
        idx = max(0, min(len(seq), 10))
    else:
        idx = first_hit_pos + 10
    return seq[:idx] + "---CUT---" + seq[idx+3: idx+3+ len(seq[idx+3:])]

def compute_delivery_success(vector: str, avg_gc: float, avg_cpg: float) -> int:
    base = {"Lipid Nanoparticles": 85, "AAV": 75, "Electroporation": 65}.get(vector, 70)
    # tweak by sequence features:
    base += -int(max(0, avg_gc-65)/5)  # too GC-rich hurts delivery
    base += -int(avg_cpg*100/10)      # high CpG hurts
    return max(5, min(98, base))

def compute_immune_risk(simulate: bool, avg_cpg: float, vector: str) -> int:
    if not simulate:
        return 0
    base = int(avg_cpg * 100 * 0.8)  # CpG heavy → immune risk
    if vector == "AAV":
        base += 5
    return max(0, min(95, base))

