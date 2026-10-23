"""
cfmd_to_csv.py: turn cFMD (curatedFoodMetagenomicData) into Neo4j import CSVs.

Run in GitHub Codespaces:
    python cfmd_to_csv.py                      # all fermented samples
    python cfmd_to_csv.py --datasets LeechJ_2020 SulaimanJ_2014   # quick test

Standard library only (no pip install needed). Output: one CSV per node type
and one per relationship type in ./import/cfmd/, plus counts.json.

Source: cFMD v1.3.2, SegataLab, CC BY 4.0
https://github.com/SegataLab/cFMD  (Carlino et al. 2024, Cell, doi 10.1016/j.cell.2024.07.039)
Every row carries a `source` column with the exact file URL it came from.
"""

import argparse
import csv
import io
import json
import re
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "v1.3.2"
BASE = f"https://raw.githubusercontent.com/SegataLab/cFMD/{VERSION}"
URL_META = f"{BASE}/cFMD_metadata.tsv"
URL_MAGS = f"{BASE}/cFMD_mags_list.tsv"
URL_DATASETS = f"{BASE}/cFMD_datasets.tsv"
CFMD_DOI = "https://doi.org/10.1016/j.cell.2024.07.039"

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus"]
PROFILE_PREFIX = {"k": "superkingdom", "p": "phylum", "c": "class",
                  "o": "order", "f": "family", "g": "genus", "s": "species", "t": "sgb"}


def profile_url(dataset):
    return f"{BASE}/cFMD_data/{dataset}/{dataset}_taxonomic_profiles.tsv"


def fetch_tsv(url):
    """Download a TSV and return a list of rows (lists). Retries 3 times."""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                text = r.read().decode("utf-8")
            return list(csv.reader(io.StringIO(text), delimiter="\t"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 * (attempt + 1))
        except urllib.error.URLError:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Could not download {url}")


def as_dicts(rows):
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:] if r]


def clean_date(s):
    """ISO date (YYYY-MM-DD) only when the raw value is unambiguous, else ''.
    cFMD mixes 2016-11-27, 11/27/2016, 2016-11, 2016 and timestamps; the raw
    string is kept in collectionDateRaw, and the year in collectionYear."""
    s = (s or "").strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}(T.*)?", s):
        return s[:10]
    return ""


def year_of(s):
    m = re.search(r"(19|20)\d{2}", s or "")
    return m.group(0) if m else ""


def to_int(s):
    try:
        return str(int(float(s)))
    except (TypeError, ValueError):
        return ""


def pretty(name):
    # cFMD writes species with underscores: Lactiplantibacillus_plantarum
    return name.replace("_", " ")


def is_placeholder(name):
    # MetaPhlAn/cFMD codes for unnamed groups: GGB123 (genus), FGB, OFGB, CFGB,
    # SGB123, EUK..., and "<group>_unclassified"
    return bool(re.match(r"(SGB|EUK|[A-Z]{0,3}GB)\d", name)) or "unclassified" in name


def known(name):
    # cFMD writes 'unknown' or leaves a rank empty when it can't assign it
    return bool(name) and name.strip().lower() not in ("unknown", "na", "none")


def taxon_id(rank, name):
    return f"{rank}:{name}"


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="import/cfmd")
    ap.add_argument("--datasets", nargs="*", help="only these datasets (test runs)")
    ap.add_argument("--include-nonfermented", action="store_true")
    ap.add_argument("--sleep", type=float, default=0.5, help="pause between downloads, s")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # 1. Metadata: which samples we keep
    meta = as_dicts(fetch_tsv(URL_META))
    keep = [m for m in meta
            if m["filtered"] == "no"
            and (args.include_nonfermented or m["fermented/non-fermented"] == "F")
            and (not args.datasets or m["dataset_name"] in args.datasets)]
    samples = {f'{m["dataset_name"]}__{m["sample_id"]}': m for m in keep}
    datasets_used = sorted({m["dataset_name"] for m in keep})
    print(f"Samples kept: {len(samples)} from {len(datasets_used)} datasets")

    # 2. Dataset table
    ds_rows = {d["Dataset"]: d for d in as_dicts(fetch_tsv(URL_DATASETS)) if d.get("Dataset")}

    # 3. MAG list
    mags_all = as_dicts(fetch_tsv(URL_MAGS))
    mags = [g for g in mags_all if f'{g["dataset_id"]}__{g["sample_id"]}' in samples]
    print(f"MAGs kept: {len(mags)}")

    # SGB info and lineage from the MAG list
    sgb = {}            # sgbId -> dict
    lineage = {}        # sgbId -> {rank: name}
    for g in mags:
        s = g["SGB_id"]
        if not s:
            continue
        sgb.setdefault(s, {"name": g["species"], "status": g["unknown"],
                           "assignedRank": g["assigned_taxonomy_level"], "source": URL_MAGS})
        lineage.setdefault(s, {r: g[r] for r in RANKS if known(g[r])})

    # 4. Taxonomic profiles, one file per dataset
    contains = []       # (sampleUid, sgbId, relAbundance, source)
    profile_lineage = {}
    missing_profiles = []
    for i, ds in enumerate(datasets_used, 1):
        url = profile_url(ds)
        rows = fetch_tsv(url)
        time.sleep(args.sleep)
        if rows is None:
            missing_profiles.append(ds)
            print(f"  [{i}/{len(datasets_used)}] {ds}: no profile file")
            continue
        header = rows[0]
        cols = []
        for c in header[1:]:
            uid = c if c in samples else f"{ds}__{c}"
            cols.append(uid if uid in samples else None)
        for r in rows[1:]:
            if not r or "|t__" not in r[0]:
                continue
            parts = dict(p.split("__", 1) for p in r[0].split("|"))
            s = parts["t"]
            profile_lineage.setdefault(s, (parts, url))
            for uid, val in zip(cols, r[1:]):
                if uid and val and float(val) > 0:
                    contains.append((uid, s, float(val), url))
        print(f"  [{i}/{len(datasets_used)}] {ds}: profile ok")

    for s, (parts, url) in profile_lineage.items():
        if s in sgb:
            continue
        lineage[s] = {rank: parts[k] for k, rank in PROFILE_PREFIX.items()
                      if rank in RANKS and known(parts.get(k))}
        sgb[s] = {"name": parts.get("s", s), "status": "", "assignedRank": "", "source": url}

    # 5. Reconcile names. cFMD mixes old and new names for the same group
    #    (Firmicutes/Bacillota, Proteobacteria/Pseudomonadota, ...), even inside
    #    the MAG list. Where one SGB carries 2 names for the same rank, treat them
    #    as synonyms; keep the name used most often in the MAG list.
    parent_of = {}                       # union-find over (rank, name)
    def find(x):
        while parent_of.get(x, x) != x:
            x = parent_of[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent_of[ra] = rb
    for s, (parts, _) in profile_lineage.items():
        mag_lin = {r: g for r, g in lineage.get(s, {}).items()}
        for k, rank in PROFILE_PREFIX.items():
            a, b = parts.get(k), mag_lin.get(rank)
            if rank in RANKS and known(a) and known(b) and a != b:
                union((rank, a), (rank, b))
    usage = Counter((r, g[r]) for g in mags for r in RANKS if known(g[r]))
    canonical = {}

    def apply_synonyms():
        groups = defaultdict(list)
        for key in set(parent_of) | set(parent_of.values()):
            groups[find(key)].append(key)
        for members in groups.values():
            best = max(members, key=lambda m: (usage.get(m, 0), m[1]))
            for m in members:
                canonical[m] = best[1]
        for lin in lineage.values():
            for rank, name in lin.items():
                lin[rank] = canonical.get((rank, name), name)

    # 6. Taxon hierarchy. A tree needs exactly 1 parent per taxon.
    def build_tree():
        taxa, votes, links = {}, defaultdict(Counter), []
        for s, lin in lineage.items():
            prev = None
            for rank in RANKS:
                name = lin.get(rank)
                if not name:
                    continue
                tid = taxon_id(rank, name)
                taxa.setdefault(tid, (rank, name, sgb[s]["source"]))
                if prev:
                    votes[tid][prev] += 1
                prev = tid
            if prev:
                links.append((s, prev, sgb[s]["source"]))
        return taxa, votes, links

    # A class, order or family listed under 2 differently named parents of the
    # same rank (Bacteroidota/Bacteroidetes) is evidence those parents are
    # synonyms: merge them, rebuild, repeat. Genus-level conflicts are NOT merged
    # (2 families can be genuinely different); they go to a majority vote.
    for _ in range(5):
        apply_synonyms()
        taxa, votes, sgb_in_taxon = build_tree()
        merged = False
        for child, c in votes.items():
            if child.split(":")[0] in ("class", "order", "family") and len(c) > 1:
                ps = [tuple(p.split(":", 1)) for p in c]
                for other in ps[1:]:
                    if find(ps[0]) != find(other):
                        union(ps[0], other)
                        merged = True
        if not merged:
            break

    taxon_parent = set()
    parent_conflicts = {}
    for child, c in votes.items():
        top = c.most_common()
        taxon_parent.add((child, top[0][0]))
        if len(top) > 1:
            parent_conflicts[child] = {p: v for p, v in top}

    # 7. Write CSVs
    n = {}
    n["datasets.csv"] = write_csv(out / "datasets.csv",
        ["datasetId", "nSamples", "accession", "cfmdVersion", "reference", "doi", "source"],
        [(d, ds_rows.get(d, {}).get("# Samples", ""), ds_rows.get(d, {}).get("Accession number", ""),
          ds_rows.get(d, {}).get("Version", ""), ds_rows.get(d, {}).get("Reference", ""), CFMD_DOI, URL_DATASETS)
         for d in datasets_used])
    n["food_types.csv"] = write_csv(out / "food_types.csv",
        ["cfmdType", "category", "fermented", "source"],
        sorted({(m["type"], m["category"], str(m["fermented/non-fermented"] == "F").lower(), URL_META)
                for m in keep}))
    n["samples.csv"] = write_csv(out / "samples.csv",
        ["sampleUid", "sampleId", "category", "cfmdType", "subtype", "commercialName", "country",
         "collectionDate", "collectionYear", "collectionDateRaw", "fermented", "studyAccession",
         "nMagsHqProk", "source"],
        [(uid, m["sample_id"], m["category"], m["type"], m["subtype"], m["commercial_name"], m["country"],
          clean_date(m["collection_date"]), year_of(m["collection_date"]), m["collection_date"],
          str(m["fermented/non-fermented"] == "F").lower(),
          m["study_accession"], to_int(m["n_MAGs_HQ_prok"]), URL_META) for uid, m in samples.items()])
    n["sgbs.csv"] = write_csv(out / "sgbs.csv",
        ["sgbId", "name", "status", "assignedRank", "isPlaceholder", "source"],
        [(s, pretty(v["name"]), v["status"], v["assignedRank"], str(is_placeholder(v["name"])).lower(), v["source"])
         for s, v in sorted(sgb.items())])
    n["taxa.csv"] = write_csv(out / "taxa.csv",
        ["taxonId", "rank", "name", "isPlaceholder", "source"],
        [(t, r, pretty(nm), str(is_placeholder(nm)).lower(), src) for t, (r, nm, src) in sorted(taxa.items())])
    n["mags.csv"] = write_csv(out / "mags.csv",
        ["magId", "completeness", "contamination", "genomeSize", "nContigs", "n50", "gcContent", "source"],
        [(g["MAG_id"], g["completeness"], g["contamination"], g["genome_size"], g["n_contigs"],
          g["N50"], g["GC_content"], URL_MAGS) for g in mags])
    # integrity check: every relationship end must exist as a node
    node_keys = {"sample": set(samples), "sgb": set(sgb), "taxon": set(taxa)}
    bad = sum(1 for c in contains if c[0] not in node_keys["sample"] or c[1] not in node_keys["sgb"])
    bad += sum(1 for s_, t, _ in sgb_in_taxon if t not in node_keys["taxon"])
    if bad:
        print(f"WARNING: {bad} relationship rows point to missing nodes")

    n["rel_sample_from_dataset.csv"] = write_csv(out / "rel_sample_from_dataset.csv",
        ["sampleUid", "datasetId", "source"],
        [(uid, m["dataset_name"], URL_META) for uid, m in samples.items()])
    n["rel_sample_of_food_type.csv"] = write_csv(out / "rel_sample_of_food_type.csv",
        ["sampleUid", "cfmdType", "source"],
        [(uid, m["type"], URL_META) for uid, m in samples.items()])
    n["rel_sample_contains_sgb.csv"] = write_csv(out / "rel_sample_contains_sgb.csv",
        ["sampleUid", "sgbId", "relAbundance", "source"], contains)
    n["rel_mag_assembled_from_sample.csv"] = write_csv(out / "rel_mag_assembled_from_sample.csv",
        ["magId", "sampleUid", "source"],
        [(g["MAG_id"], f'{g["dataset_id"]}__{g["sample_id"]}', URL_MAGS) for g in mags])
    n["rel_mag_assigned_to_sgb.csv"] = write_csv(out / "rel_mag_assigned_to_sgb.csv",
        ["magId", "sgbId", "source"],
        [(g["MAG_id"], g["SGB_id"], URL_MAGS) for g in mags if g["SGB_id"]])
    n["rel_sgb_in_taxon.csv"] = write_csv(out / "rel_sgb_in_taxon.csv",
        ["sgbId", "taxonId", "source"], sgb_in_taxon)
    n["rel_taxon_part_of.csv"] = write_csv(out / "rel_taxon_part_of.csv",
        ["childTaxonId", "parentTaxonId", "source"],
        [(c, p, taxa[c][2]) for c, p in sorted(taxon_parent)])

    nodes = sum(v for k, v in n.items() if not k.startswith("rel_"))
    rels = sum(v for k, v in n.items() if k.startswith("rel_"))
    report = {"cfmdVersion": VERSION, "rows_per_file": n, "nodes_total": nodes, "rels_total": rels,
              "missing_profiles": missing_profiles,
              "synonyms_merged": {f"{r}:{a}": b for (r, a), b in sorted(canonical.items()) if a != b},
              "taxa_with_conflicting_parents": parent_conflicts}
    (out / "counts.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
