#!/usr/bin/env python3
"""Build a LINCS L1000 signature matrix (compounds x genes) for the perturbational
stage, from the public Enrichr LINCS L1000 chemical-perturbation gene-set libraries.

No login required. Two libraries are used:
  LINCS_L1000_Chem_Pert_up    (genes up after each perturbation)
  LINCS_L1000_Chem_Pert_down  (genes down after each perturbation)

Each perturbation term becomes a signed signature: +1 for up genes, -1 for down.
Signatures are aggregated per compound (mean over all its signatures), which is
exactly the cross-signature recurrence the framework rewards.

Outputs:
  <out>/lincs_signatures.csv   compounds x genes signed matrix
  <out>/drug_meta.csv          compound,target_class,risk_class (filled where known)

Usage:  python scripts/build_lincs.py --out data/lincs
"""
from __future__ import annotations
import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="

# Drug-NAMED sources (terms begin with the perturbagen name). Tried in order.
SOURCES = {
    "geo": {"up": ["Drug_Perturbations_from_GEO_up"],
            "down": ["Drug_Perturbations_from_GEO_down"]},
    "l1000": {"up": ["LINCS_L1000_Chem_Pert_up"],
              "down": ["LINCS_L1000_Chem_Pert_down"]},
}

# plate / experiment codes that are NOT drug names (L1000 keying artefacts)
PLATE_RE = re.compile(r"^(brd|ljp|cpc|rep|dos|trcn|aby|kd|oe)[-_]?\d", re.I)


def fetch_gmt(candidates: list[str]) -> str:
    last = None
    for name in candidates:
        try:
            print(f"  downloading {name} ...")
            r = requests.get(BASE + name, timeout=300)
            r.raise_for_status()
            if r.text.strip():
                return r.text
        except Exception as e:  # noqa
            last = e
            print(f"    failed: {e}")
    raise SystemExit(f"Could not download any of {candidates}: {last}")


_STOP = re.compile(r"^(db\d+|gse\d+|gds\d+|\d+.*|human|mouse|rat|sample|cid|sid|"
                   r"na|compound|drug|up|down|dn|\d+x.*)$", re.I)


def compound_of(term: str) -> str:
    """Extract the clean drug name from a GEO drug-perturbation term such as
    'hydrocortisone db00741 human gse7890 sample 2751' or
    'valproic acid db00313 mouse gse9247 sample 2838'. Keep leading word tokens
    until an identifier token (DBID, number, species, 'sample', 'gse...') appears,
    so every sample of one drug collapses to the same name (aggregated later)."""
    name = []
    for tok in re.split(r"[\s,]+", term.strip()):
        if not tok:
            continue
        if _STOP.match(tok) or re.search(r"\d", tok):
            break
        name.append(tok)
    return " ".join(name).lower().strip()


def is_drug_name(name: str) -> bool:
    if not name or len(name) < 3:
        return False
    if PLATE_RE.match(name):
        return False
    if not re.search(r"[aeiou]", name):   # codes rarely have vowels
        return False
    return bool(re.match(r"^[a-z]", name))


def parse_gmt(text: str) -> dict[str, set]:
    sig = {}
    for line in text.strip().splitlines():
        parts = line.split("\t")
        term = parts[0]
        genes = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
        if genes:
            sig[term] = genes
    return sig


def build_matrix(up: dict, down: dict, max_genes: int = 4000):
    # per-compound aggregation
    comp_up: dict[str, list[set]] = defaultdict(list)
    comp_dn: dict[str, list[set]] = defaultdict(list)
    for term, genes in up.items():
        c = compound_of(term)
        if is_drug_name(c):
            comp_up[c].append(genes)
    for term, genes in down.items():
        c = compound_of(term)
        if is_drug_name(c):
            comp_dn[c].append(genes)

    compounds = sorted(set(comp_up) | set(comp_dn))
    # restrict to the most frequently perturbed genes to keep the matrix compact
    gene_freq: dict[str, int] = defaultdict(int)
    for genes_list in list(comp_up.values()) + list(comp_dn.values()):
        for genes in genes_list:
            for g in genes:
                gene_freq[g] += 1
    top_genes = [g for g, _ in sorted(gene_freq.items(),
                                      key=lambda kv: kv[1], reverse=True)[:max_genes]]
    gidx = {g: i for i, g in enumerate(top_genes)}

    mat = np.zeros((len(compounds), len(top_genes)), dtype=float)
    for r, c in enumerate(compounds):
        n = 0
        for genes in comp_up.get(c, []):
            for g in genes:
                if g in gidx:
                    mat[r, gidx[g]] += 1.0
            n += 1
        for genes in comp_dn.get(c, []):
            for g in genes:
                if g in gidx:
                    mat[r, gidx[g]] -= 1.0
            n += 1
        if n:
            mat[r] /= n  # mean signed recurrence
    return pd.DataFrame(mat, index=compounds, columns=top_genes)


def make_drug_meta(compounds):
    """Fill target_class / risk_class from config/drug_safety.yaml where the
    compound name matches a known entry; otherwise leave blank."""
    import yaml
    safety = yaml.safe_load(open(ROOT / "config" / "drug_safety.yaml"))
    # a light lookup table for the demo panel + common repositioning drugs
    known = {
        "nintedanib": ("pdgf", ""), "pirfenidone": ("tgfb", ""),
        "saracatinib": ("src", ""), "imatinib": ("pdgf", ""),
        "dasatinib": ("src", "cytotoxic_chemo"), "gefitinib": ("pdgf", ""),
        "lapatinib": ("pdgf", ""), "foretinib": ("pdgf", ""), "pd173074": ("pdgf", ""),
        "sirolimus": ("mtor", ""), "rapamycin": ("mtor", ""),
        "vorinostat": ("hdac", ""), "valproic acid": ("hdac", ""),
        "ruxolitinib": ("jak", "immunosuppressant_high"),
        "fasudil": ("rho_rock", ""), "losartan": ("tgfb", ""),
        # interpretable repurposing classes (fibrosis-relevant targets)
        "rosiglitazone": ("pparg", ""), "pioglitazone": ("pparg", ""),
        "troglitazone": ("pparg", ""), "metformin": ("ampk", ""),
        "atorvastatin": ("hmgcr", ""), "simvastatin": ("hmgcr", ""),
        "quercetin": ("senolytic", ""), "resveratrol": ("sirtuin", ""),
        "tretinoin": ("retinoid", ""), "bexarotene": ("retinoid", ""),
        "curcumin": ("nfkb", ""),
        # deprioritized: environmental / genotoxic / cytotoxic-chemo / steroid
        "haloperidol": ("", "neuropsychiatric"), "clozapine": ("", "neuropsychiatric"),
        "thalidomide": ("", "neuropsychiatric"), "colchicine": ("", "cytotoxic_chemo"),
    }
    CHEMO = {"cisplatin", "carboplatin", "doxorubicin", "paclitaxel", "irinotecan",
             "methotrexate", "bleomycin", "bortezomib", "decitabine", "puromycin",
             "streptozocin", "mycophenolic acid", "vemurafenib", "plx4032"}
    ENV = {"formaldehyde", "bpde", "cadmium", "cadmium dichloride", "trichloroethylene",
           "sodium arsenite", "sodium dichromate", "nickel chloride", "cobalt dichloride",
           "cobalt dichloride hexahydrate", "bisphenol a", "perfluorooctanoic acid",
           "mehp", "pbde", "cyfluthrin", "chlorpyrifos", "deltamethrin", "permethrin",
           "monuron", "titanium dioxide", "hydrogen peroxide", "cigarette smoke",
           "ethanol", "heroin", "methamphetamine", "morphine", "nicotine",
           "trichloroethylene", "diethylstilbestrol"}
    GENOTOX = {"bpde", "sodium arsenite", "diethylstilbestrol", "cadmium"}
    HORMONE = {"estradiol", "ethinyl estradiol", "testosterone", "r1881", "genistein",
               "letrozole", "anastrozole", "tibolone", "abiraterone"}
    rows = []
    for c in compounds:
        if c in known:
            t, r = known[c]
        elif c in GENOTOX:
            t, r = "", "genotoxic"
        elif c in ENV:
            t, r = "", "environmental"
        elif c in CHEMO:
            t, r = "", "cytotoxic_chemo"
        elif c in HORMONE:
            t, r = "", "hormone"
        else:
            t, r = "", ""
        rows.append({"compound": c, "target_class": t, "risk_class": r})
    return pd.DataFrame(rows).set_index("compound")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/lincs")
    ap.add_argument("--max-genes", type=int, default=4000)
    ap.add_argument("--source", choices=list(SOURCES), default="geo",
                    help="'geo' = drug-NAMED perturbations (recommended); "
                         "'l1000' = L1000 chem-pert (often plate-coded IDs)")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    libs = SOURCES[args.source]
    print(f"Building drug signature matrix from Enrichr source '{args.source}'")
    up = parse_gmt(fetch_gmt(libs["up"]))
    dn = parse_gmt(fetch_gmt(libs["down"]))
    print(f"  {len(up)} up-terms, {len(dn)} down-terms")
    mat = build_matrix(up, dn, max_genes=args.max_genes)
    mat.to_csv(out / "lincs_signatures.csv")
    meta = make_drug_meta(mat.index)
    meta.to_csv(out / "drug_meta.csv")
    print(f"  wrote {out/'lincs_signatures.csv'}  "
          f"({mat.shape[0]} compounds x {mat.shape[1]} genes)")
    print(f"  wrote {out/'drug_meta.csv'}")


if __name__ == "__main__":
    main()
