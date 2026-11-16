"""Systematic verification of CIRCE's Morillo 2025 male neuropeptide layer against Haque et al. 2025.

Two parts, both over the neuropeptide ligands (FLP/NLP/INS) and NPP-GPCR receptors of the connectome:

  PART A -- concordance (does Haque corroborate the Morillo calls CIRCE ingested?)
    A1. AUC: how well Haque's continuous per-class mean expression (Tables S7 herm / S8 male)
        separates the neuron classes Morillo called positive from those it called negative.
    A2. Cohen's kappa + agreement between Morillo's binary calls and Haque's binary calls (Table S6).
    Headline over Haque's 62 high-confidence sex-shared clusters (Table S5); all shared classes as a
    sensitivity check. Ligands and receptors reported separately (receptors are the connectome's
    weak, edge-bearing layer).

  PART B -- remodeling cross-check (the sex comparison Morillo's depth confound forced us to drop)
    B1. Haque's OWN male-vs-herm neuropeptide bias, platform-matched, from S7/S8 means and the S3
        differential-expression table -- tests the paper's claim that most neuropeptides are
        male-biased, independently of Morillo.
    B2. Per-gene DIRECTION agreement between Morillo (confounded) and Haque (clean): does the Morillo
        depth confound flip individual genes' sex-bias direction, or only inflate herm globally?

Haque supplementary tables are large untrusted third-party downloads and are NOT vendored. Give their
directory via $HAQUE_DIR or argv[1]; the committed results.json captures the numbers. Deterministic
(sorted; no RNG). Run from a safe cwd:

    HAQUE_DIR=/path/to/PMID40644535 uv run python -I analysis/haque_morillo_validation/validate.py
"""

import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import openpyxl

REPO = Path(__file__).resolve().parents[2]
HAQUE = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("HAQUE_DIR", ""))
MALE = REPO / "data" / "morillo-2025-male"
OUT = Path(__file__).parent / "results.json"

ALIAS = {"ctdh-1": "nlp-73", "srlf-1": "nlp-77"}  # Haque canonical -> Morillo symbol


def gene_set(path):
    return {r["symbol"]: r["wbgene"] for r in csv.DictReader(open(path))}


LIGANDS = gene_set(MALE / "morillo_genes.csv")
RECEPTORS = gene_set(MALE / "morillo_receptor_genes.csv")
TARGETS = set(LIGANDS) | set(RECEPTORS) | set(ALIAS)

# cell/label -> CIRCE class (Haque's column namespace is class-level)
conn = json.load(open(REPO / "outputs" / "connectome.json"))
CLASS = {c["name"]: (c.get("cell_class") or c["name"]) for c in conn["cells"]}


def load_morillo(path):
    out = defaultdict(set)
    for r in csv.DictReader(open(path)):
        out[r["gene_symbol"]].add(CLASS.get(r["cell"], r["cell"]))
    return out


MOR = {
    ("hermaphrodite", "lig"): load_morillo(MALE / "morillo_herm_expression.csv"),
    ("male", "lig"): load_morillo(MALE / "morillo_male_expression.csv"),
    ("hermaphrodite", "rec"): load_morillo(MALE / "morillo_herm_receptor_expression.csv"),
    ("male", "rec"): load_morillo(MALE / "morillo_male_receptor_expression.csv"),
}


def _sheet(path, name):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[name]
    it = ws.iter_rows(values_only=True)
    header = [str(c).strip() if c is not None else "" for c in next(it)]
    for row in it:
        yield header, row
    wb.close()


def load_mean(path, sheet):
    """gene -> {class: mean expression} for TARGET genes only."""
    data, cols = {}, None
    for header, row in _sheet(path, sheet):
        cols = header[1:]
        g = row[0]
        if g is None or str(g).strip() not in TARGETS:
            continue
        data[str(g).strip()] = {
            cols[i]: float(v) for i, v in enumerate(row[1:]) if isinstance(v, (int, float))
        }
    return set(cols or []), data


def load_binary(path, sheet):
    """gene -> {column: 1/0} for TARGET genes (columns are e.g. ADA_Herm / ADA_Male)."""
    data, cols = {}, None
    for header, row in _sheet(path, sheet):
        cols = header[1:]
        g = row[0]
        if g is None or str(g).strip() not in TARGETS:
            continue
        data[str(g).strip()] = {
            cols[i]: (1 if (isinstance(v, (int, float)) and v > 0) else 0)
            for i, v in enumerate(row[1:])
        }
    return data


HC_H, MEAN_H = load_mean(HAQUE / "adv9106_table_s7.xlsx", "Herm.Exp.Mean")
HC_M, MEAN_M = load_mean(HAQUE / "adv9106_table_s8.xlsx", "Male.Exp.Mean")
BIN = load_binary(HAQUE / "adv9106_table_s6.xlsx", "Sheet 1 - 240215_threshold_by_t")
MEAN = {"hermaphrodite": (HC_H, MEAN_H), "male": (HC_M, MEAN_M)}

# 62 high-confidence sex-shared clusters (Table S5)
HICONF = set()
for _h, row in _sheet(HAQUE / "adv9106_table_s5.xlsx", "Sheet1"):
    if row[0]:
        HICONF.add(str(row[0]).strip())


def hexp(gene, vec_by_gene):
    return vec_by_gene.get(gene) or vec_by_gene.get(ALIAS.get(gene, gene)) or {}


def auc(pos, neg):
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def kappa(pairs):
    """Cohen's kappa over [(a,b)] of 0/1 labels."""
    n = len(pairs)
    if not n:
        return None
    po = sum(a == b for a, b in pairs) / n
    pa1 = sum(a for a, _ in pairs) / n
    pb1 = sum(b for _, b in pairs) / n
    pe = pa1 * pb1 + (1 - pa1) * (1 - pb1)
    return (po - pe) / (1 - pe) if pe != 1 else 1.0


def concordance(classes, restrict_label):
    res = {}
    for sex in ("male", "hermaphrodite"):
        hcols, hmean = MEAN[sex]
        usable = classes & hcols
        for kind, genes in (("lig", LIGANDS), ("rec", RECEPTORS)):
            mcalls = MOR[(sex, kind)]
            pos, neg, per_gene, kpairs, discordant = [], [], [], [], []
            for g in sorted(genes):
                vec = {k: v for k, v in hexp(g, hmean).items() if k in usable}
                if not vec:
                    continue
                poslab = mcalls.get(g, set())
                gp = [v for k, v in vec.items() if k in poslab]
                gn = [v for k, v in vec.items() if k not in poslab]
                a = auc(gp, gn)
                if a is not None:
                    per_gene.append((g, a))
                    pos += gp
                    neg += gn
                    if a < 0.5:
                        discordant.append((g, round(a, 3)))
                # binary kappa vs Haque S6
                suffix = "_Male" if sex == "male" else "_Herm"
                hb = hexp(g, BIN)
                for cl in usable:
                    col = f"{cl}{suffix}"
                    if col in hb:
                        kpairs.append((1 if cl in poslab else 0, hb[col]))
            res[f"{sex}.{kind}"] = {
                "genes_compared": len(per_gene),
                "classes_used": len(usable),
                "pooled_auc": round(auc(pos, neg), 3) if pos and neg else None,
                "mean_per_gene_auc": (
                    round(sum(a for _, a in per_gene) / len(per_gene), 3) if per_gene else None
                ),
                "binary_kappa": round(kappa(kpairs), 3) if kpairs else None,
                "binary_agreement": (
                    round(sum(a == b for a, b in kpairs) / len(kpairs), 3) if kpairs else None
                ),
                "discordant_genes_auc_below_0.5": sorted(discordant, key=lambda x: x[1]),
            }
    return {"restrict": restrict_label, "metrics": res}


# ---- PART A: concordance (62 high-confidence clusters; all-shared as sensitivity) ----
all_classes = HC_H | HC_M
partA = {
    "high_confidence_62": concordance(HICONF, "Haque Table S5 (62 high-confidence sex-shared)"),
    "all_shared_sensitivity": concordance(all_classes, "all annotated shared classes"),
}

# ---- PART B: remodeling cross-check ----
# B1: Haque's own male-vs-herm neuropeptide bias over the 62 clusters (platform-matched).
shared62 = HICONF & HC_H & HC_M


def haque_bias(genes):
    male_biased = herm_biased = flat = 0
    details = {}
    for g in sorted(genes):
        hm, hh = hexp(g, MEAN_M), hexp(g, MEAN_H)
        clusters = [c for c in shared62 if c in hm and c in hh]
        if not clusters:
            continue
        dm = sum(hm[c] for c in clusters) / len(clusters)
        dh = sum(hh[c] for c in clusters) / len(clusters)
        direction = "male" if dm > dh else ("herm" if dh > dm else "flat")
        details[g] = {"male_mean": round(dm, 3), "herm_mean": round(dh, 3), "dir": direction}
        male_biased += direction == "male"
        herm_biased += direction == "herm"
        flat += direction == "flat"
    return male_biased, herm_biased, flat, details


# S3 significant DE events for neuropeptide genes (male-biased = positive logFC).
s3_male = s3_herm = 0
for _h, row in _sheet(HAQUE / "adv9106_table_s3.xlsx", "Table S3"):
    g = row[0]
    if g is None or str(g).strip() not in set(LIGANDS):
        continue
    try:
        padj, lfc = float(row[6]), float(row[7])
    except (TypeError, ValueError):
        continue
    if padj < 0.05:
        s3_male += lfc > 0
        s3_herm += lfc < 0

lig_mb, lig_hb, lig_flat, lig_detail = haque_bias(LIGANDS)

# B2: per-gene direction agreement, Morillo (confounded) vs Haque (clean).
agree = disagree = 0
for g in sorted(LIGANDS):
    if g not in lig_detail:
        continue
    m_male = len([c for c in MOR[("male", "lig")].get(g, set()) if c in shared62])
    m_herm = len([c for c in MOR[("hermaphrodite", "lig")].get(g, set()) if c in shared62])
    if m_male == m_herm:
        continue  # Morillo has no directional call for this gene on the shared set
    m_dir = "male" if m_male > m_herm else "herm"
    if lig_detail[g]["dir"] in ("male", "herm"):
        agree += m_dir == lig_detail[g]["dir"]
        disagree += m_dir != lig_detail[g]["dir"]

partB = {
    "haque_ligand_sex_bias_62clusters": {
        "male_biased_genes": lig_mb,
        "herm_biased_genes": lig_hb,
        "flat": lig_flat,
        "note": "tests Haque's claim that most neuropeptides are male-biased (platform-matched)",
    },
    "haque_S3_significant_DE_events_ligands": {
        "male_biased_events_padj<0.05": s3_male,
        "herm_biased_events_padj<0.05": s3_herm,
    },
    "morillo_vs_haque_ligand_direction": {
        "agree": agree,
        "disagree": disagree,
        "agreement_rate": round(agree / (agree + disagree), 3) if (agree + disagree) else None,
        "note": "does Morillo's depth confound flip per-gene sex-bias direction, or only inflate herm globally?",
    },
}

result = {
    "haque_source": "Haque et al. 2025 Sci Adv (adv9106) supplementary tables S3/S5/S6/S7/S8",
    "part_a_concordance": partA,
    "part_b_remodeling": partB,
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
