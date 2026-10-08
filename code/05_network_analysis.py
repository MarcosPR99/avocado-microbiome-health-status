#!/usr/bin/env python3
"""
Cross-kingdom co-association network analysis (SparCC + NetworkX).

Six networks (HR, DR, HE, DE, HP, DP) were inferred from combined October
16S rRNA gene and ITS abundance tables, at ASV and genus resolution.

    ASV level:   top 500 bacterial + top 500 fungal ASVs per condition
                 (total abundance across shared samples)
    Genus level: ASVs collapsed to genus; genera retained if present in
                 >= 3 samples and with >= 50 total reads
    SparCC:      -i 20; two-sided pseudo-P-values from 99 permutations
    Edges:       |rho| > 0.70 and P <= 0.01
    Metrics:     nodes, edges, % positive/negative edges, modularity
                 (greedy modularity, unweighted), average degree, diameter
                 and average path length (largest connected component),
                 % bacteria-bacteria / bacteria-fungi / fungi-fungi edges

Bacterial and fungal features are prefixed B_ and F_. Edge tables can be
imported directly into Gephi (Source, Target, Weight).

Requirements: Python 3.8, pandas, NetworkX 3.1, and a local SparCC
installation (SparCC.py, MakeBootstraps.py, PseudoPvals.py).

Usage:
    python 05_network_analysis.py prepare --level ASV   --in16 DIR --inits DIR --out DIR
    python 05_network_analysis.py prepare --level GENUS --in16 DIR --inits DIR --out DIR
    python 05_network_analysis.py sparcc  --sparcc-dir DIR --inputs DIR --out DIR
    python 05_network_analysis.py edges   --runs DIR --out DIR --level ASV
    python 05_network_analysis.py metrics --edges DIR --out FILE --level ASV
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
from networkx.algorithms.community import greedy_modularity_communities, modularity

CONDITIONS = ["HR", "DR", "HE", "DE", "HP", "DP"]
TOP_N = 500
MIN_PREVALENCE = 3
MIN_READS = 50
ITERATIONS = 20
PERMUTATIONS = 99
RHO = 0.70
P_MAX = 0.01

NON_SAMPLE_COLUMNS = {
    "taxonomy", "Taxonomy", "Taxon", "taxon", "Confidence", "confidence",
    "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
}
UNCLASSIFIED = {"", "nan", "none", "na", "unassigned", "uncultured", "unidentified"}


# --------------------------------------------------------------------------
# Input tables
# --------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    return df.set_index(df.columns[0])


def sample_columns(df: pd.DataFrame) -> list:
    return [c for c in df.columns
            if c not in NON_SAMPLE_COLUMNS and pd.to_numeric(df[c], errors="coerce").notna().any()]


def counts(df: pd.DataFrame, samples: list) -> pd.DataFrame:
    return df[samples].apply(pd.to_numeric, errors="coerce").fillna(0)


def genus_labels(df: pd.DataFrame) -> pd.Series:
    """Genus from the taxonomy string; unclassified features keep their own ID."""
    col = next(c for c in ("Taxon", "taxon", "taxonomy", "Taxonomy") if c in df.columns)
    labels = []
    for feature_id, tax in df[col].items():
        ranks = [r.strip() for r in str(tax).split(";")]
        genus = next((r[3:] for r in reversed(ranks) if r.lower().startswith("g__")), "")
        genus = genus.strip().replace(" ", "_")
        labels.append(f"Unclassified_{feature_id}" if genus.lower() in UNCLASSIFIED else genus)
    return pd.Series(labels, index=df.index)


def prepare(args: argparse.Namespace) -> None:
    args.out.mkdir(parents=True, exist_ok=True)
    for cond in args.conditions:
        b = read_table(args.in16 / f"ASV_table_{cond}_Oct_withTax_ASV.tsv")
        f = read_table(args.inits / f"ITS_table_{cond}_Oct_withTax_ASV.tsv")
        shared = [s for s in sample_columns(b) if s in set(sample_columns(f))]

        tables = []
        for df, prefix in ((b, "B_"), (f, "F_")):
            x = counts(df, shared)
            if args.level == "ASV":
                x = x.loc[x.sum(axis=1).sort_values(ascending=False).index[:TOP_N]]
            else:
                x = x.groupby(genus_labels(df)).sum()
            x.index = [prefix + str(i) for i in x.index]
            tables.append(x)

        table = pd.concat(tables)
        if args.level == "GENUS":
            table = table[((table > 0).sum(axis=1) >= MIN_PREVALENCE) & (table.sum(axis=1) >= MIN_READS)]

        out = args.out / f"{cond}_{args.level}.tsv"
        table.to_csv(out, sep="\t", index_label="FeatureID")
        print(f"{cond}: {table.shape[0]} features, {len(shared)} samples -> {out}")


# --------------------------------------------------------------------------
# SparCC
# --------------------------------------------------------------------------

def run(script: Path, *params: str) -> None:
    subprocess.run([sys.executable, str(script), *params], check=True)


def sparcc(args: argparse.Namespace) -> None:
    for table in sorted(args.inputs.glob("*.tsv")):
        run_dir = args.out / table.stem
        perm_dir, perm_cor_dir = run_dir / "perm", run_dir / "perm_cor"
        perm_cor_dir.mkdir(parents=True, exist_ok=True)
        perm_dir.mkdir(parents=True, exist_ok=True)

        run(args.sparcc_dir / "SparCC.py", str(table), "-i", str(ITERATIONS),
            f"--cor_file={run_dir / 'cor_sparcc.out'}")
        run(args.sparcc_dir / "MakeBootstraps.py", str(table), "-n", str(PERMUTATIONS),
            "-t", "permutation_#.txt", "-p", f"{perm_dir}/")
        for i in range(PERMUTATIONS):
            run(args.sparcc_dir / "SparCC.py", str(perm_dir / f"permutation_{i}.txt"),
                "-i", str(ITERATIONS), f"--cor_file={perm_cor_dir / f'permutation_{i}_cor.out'}")
        run(args.sparcc_dir / "PseudoPvals.py", str(run_dir / "cor_sparcc.out"),
            str(perm_cor_dir / "permutation_#_cor.out"), str(PERMUTATIONS),
            "-o", str(run_dir / "pvals_two_sided.txt"), "-t", "two_sided")
        print(f"{table.stem}: SparCC done")


# --------------------------------------------------------------------------
# Edges and network metrics
# --------------------------------------------------------------------------

def edges(args: argparse.Namespace) -> None:
    args.out.mkdir(parents=True, exist_ok=True)
    for cond in args.conditions:
        run_dir = args.runs / f"{cond}_{args.level}"
        rho = pd.read_csv(run_dir / "cor_sparcc.out", sep="\t", index_col=0)
        pval = pd.read_csv(run_dir / "pvals_two_sided.txt", sep="\t", index_col=0)

        taxa = list(rho.columns)
        rows = []
        for i, a in enumerate(taxa):
            for b in taxa[i + 1:]:
                r, p = float(rho.at[a, b]), float(pval.at[a, b])
                if abs(r) > RHO and p <= P_MAX:
                    rows.append((a, b, r, p))

        out = args.out / f"{cond}_{args.level}_edges.tsv"
        pd.DataFrame(rows, columns=["Source", "Target", "Weight", "Pvalue"]).to_csv(out, sep="\t", index=False)
        print(f"{cond}: {len(rows)} edges -> {out}")


def network_metrics(df: pd.DataFrame) -> dict:
    G = nx.from_pandas_edgelist(df, "Source", "Target")
    n_nodes, n_edges = G.number_of_nodes(), G.number_of_edges()
    if n_edges == 0:
        return {"Nodes": n_nodes, "Edges": 0}

    lcc = G.subgraph(max(nx.connected_components(G), key=len))
    kind = df["Source"].str[0] + df["Target"].str[0]

    return {
        "Nodes": n_nodes,
        "Edges": n_edges,
        "Positive_%": 100 * (df["Weight"] > 0).mean(),
        "Negative_%": 100 * (df["Weight"] < 0).mean(),
        "Modularity": modularity(G, greedy_modularity_communities(G, weight=None), weight=None),
        "Average_degree": 2 * n_edges / n_nodes,
        "Diameter_LCC": nx.diameter(lcc),
        "Average_path_length_LCC": nx.average_shortest_path_length(lcc),
        "BB_%": 100 * (kind == "BB").mean(),
        "BF_%": 100 * kind.isin(["BF", "FB"]).mean(),
        "FF_%": 100 * (kind == "FF").mean(),
    }


def metrics(args: argparse.Namespace) -> None:
    rows = []
    for cond in args.conditions:
        df = pd.read_csv(args.edges / f"{cond}_{args.level}_edges.tsv", sep="\t")
        rows.append({"Network": cond, "Level": args.level, **network_metrics(df)})
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Metrics -> {args.out}")


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare", help="Build cross-kingdom SparCC input tables")
    p.add_argument("--level", choices=["ASV", "GENUS"], required=True)
    p.add_argument("--in16", type=Path, required=True, help="Directory with 16S ASV tables")
    p.add_argument("--inits", type=Path, required=True, help="Directory with ITS ASV tables")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=prepare)

    p = sub.add_parser("sparcc", help="Run SparCC, permutations and pseudo-P-values")
    p.add_argument("--sparcc-dir", type=Path, required=True)
    p.add_argument("--inputs", type=Path, required=True, help="Directory with tables from 'prepare'")
    p.add_argument("--out", type=Path, required=True)
    p.set_defaults(func=sparcc)

    p = sub.add_parser("edges", help="Retain edges with |rho| > 0.70 and P <= 0.01")
    p.add_argument("--runs", type=Path, required=True, help="Output directory of 'sparcc'")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--level", choices=["ASV", "GENUS"], required=True)
    p.set_defaults(func=edges)

    p = sub.add_parser("metrics", help="Compute network properties (Tables S9-S10)")
    p.add_argument("--edges", type=Path, required=True, help="Output directory of 'edges'")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--level", choices=["ASV", "GENUS"], required=True)
    p.set_defaults(func=metrics)

    for p in sub.choices.values():
        p.add_argument("--conditions", nargs="+", default=CONDITIONS)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
