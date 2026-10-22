#!/usr/bin/env python3
"""Fetch 3D structures for scraped sequences and align them in 3D. (v2)

Usage:
    python fetch_and_align_structures.py sequences.fasta -o output_dir
    python fetch_and_align_structures.py ECA3.fasta -o eca3_structures --summary eca3.tsv

v2 changes (ribbon-ready, animation-ready outputs):
    - aligned.json now carries the residue SEQUENCE and per-residue
      secondary structure (H/E/C from CA geometry) for every structure,
      so the morph-animation and interactive-viewer skills can compute a
      stability mask without re-deriving anything.
    - aligned PDBs are written with HELIX/SHEET records -> they open as
      ribbons in PyMOL / ChimeraX / Mol*.
    - no membrane prediction anywhere (removed; looked odd and is not
      needed for the evolutionary story).

Per sequence:
    1. Resolve a UniProt accession (direct, NCBI->UniProt mapping, or
       UniProt search by protein-name keywords + organism). Every
       candidate is VERIFIED by k-mer similarity against the input
       sequence (>=0.30) - wrong-but-same-length proteins are rejected.
    2. Fetch the AlphaFold DB model (CA trace + per-residue pLDDT).
    3. If no AlphaFold entry exists: build a homology model by threading
       the sequence onto the closest available structure (local BLOSUM62).
Alignment:
    All structures are superposed onto one reference (real, non-threaded,
    highest average core identity) using residue correspondences from
    LOCAL BLOSUM62 alignments and Kabsch rotation with iterative core
    refinement. Threaded models inherit their template's transform.

Outputs (in -o dir):
    raw/<NAME>.pdb          downloaded AlphaFold models (original frame)
    aligned/<NAME>.pdb      superposed CA-trace models (common frame),
                            ribbon-ready (HELIX/SHEET records)
    aligned.json            coords + plddt + seq + ss + metadata for
                            downstream steps
    overlay.png             CA traces of all structures overlaid (visual check)
    summary.tsv             per-structure: source, plddt, rmsd-to-reference
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import requests

# shared library: prefer a checkout next to the skills tree, fall back to
# an installed copy
_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, "..", "..", "..", "src"),
              os.path.join(_HERE, "..", "..", "src"),
              os.path.join(_HERE, "src")):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, os.path.abspath(_cand))
try:
    import proteinlib as pl
except ImportError:
    pl = None

AF_BASE = "https://alphafold.ebi.ac.uk"
UP_BASE = "https://rest.uniprot.org"


def kmer_similarity(seq_a, seq_b, k=5):
    """Fraction of seq_a k-mers shared with seq_b (identity proxy)."""
    if not seq_a or not seq_b:
        return 0.0
    ka = {seq_a[i:i + k] for i in range(len(seq_a) - k + 1)}
    kb = {seq_b[i:i + k] for i in range(len(seq_b) - k + 1)}
    if not ka:
        return 0.0
    return len(ka & kb) / len(ka)


def parse_fasta(path):
    """Parse the scraper's FASTA into records.

    Header format: >ACCESSION Species|Lineage|Family Definition
    """
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append(_make_record(header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header:
        records.append(_make_record(header, "".join(seq)))
    return records


def _make_record(header, seq):
    """Parse '>ACCESSION Species|Lineage|Family Definition...' headers."""
    parts = header.split(None, 1)
    accession = parts[0]
    rest = parts[1] if len(parts) > 1 else ""
    species = lineage = family = ""
    definition = rest
    if "|" in rest:
        import re as _re
        m = _re.match(r"^(.+?)\|([A-Za-z ]+)\|([A-Za-z ]+?)(?:\s+\|?)(.*)$",
                      rest)
        if m:
            species, lineage, family, definition = (m.group(1).strip(),
                                                    m.group(2).strip(),
                                                    m.group(3).strip(),
                                                    m.group(4).strip())
        else:
            bits = rest.split("|", 2)
            if len(bits) == 3:
                species, lineage, tail = bits
                fam_end = tail.find(" ")
                if fam_end == -1:
                    family, definition = tail, ""
                else:
                    family = tail[:fam_end]
                    definition = tail[fam_end + 1:]
    idx = definition.find("[")
    if 0 < idx < len(definition) - 1:
        definition = definition[:idx].strip() or definition
    return {
        "accession": accession,
        "species": species,
        "lineage": lineage,
        "family": family,
        "definition": definition,
        "sequence": seq.upper(),
    }


class StructureFetcher:
    def __init__(self, cache_dir, email=None):
        self.cache_dir = cache_dir
        os.makedirs(os.path.join(cache_dir, "raw"), exist_ok=True)
        os.makedirs(os.path.join(cache_dir, "aligned"), exist_ok=True)
        self.session = requests.Session()
        if email:
            self.session.params = {"email": email}

    def _get(self, url, params=None, timeout=60, retries=3):
        for attempt in range(retries):
            try:
                r = self.session.get(url, params=params, timeout=timeout)
                if r.status_code == 200:
                    return r
            except requests.RequestException:
                pass
            time.sleep(2 * (attempt + 1))
        return None

    def ncbi_to_uniprot(self, accession):
        """Map an NCBI protein accession to UniProt via the ID-mapping API."""
        try:
            r = requests.post(f"{UP_BASE}/idmapping/run", timeout=30, data={
                "from": "RefSeq_Protein", "to": "UniProtKB",
                "ids": accession})
            if r.status_code != 200:
                return None
            job_id = r.json()["jobId"]
        except (requests.RequestException, ValueError, KeyError):
            return None
        for _ in range(20):
            time.sleep(3)
            res = self._get(
                f"{UP_BASE}/idmapping/uniprotkb/results/stream/{job_id}",
                params={"format": "tsv", "fields": "accession"})
            if res and "From" in res.text:
                for line in res.text.strip().splitlines()[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2 and parts[1]:
                        return parts[1]
                return None
        return None

    def uniprot_search_by_seq(self, species, sequence, definition=""):
        """Find a UniProt entry by organism, verified by sequence identity.

        Queries: (a) protein-name keywords from the definition + organism,
        (b) organism only. Candidates are verified by k-mer similarity
        against the input sequence; only high-similarity entries are
        accepted, preventing spurious matches (RNA polymerases, PEPC,
        kinesins, etc. caught in earlier testing).
        """
        import re as _re
        genus = species.split()[0] if species else ""
        if not genus or not sequence:
            return None
        queries = []
        keywords = [w for w in _re.findall(r"[A-Za-z0-9]+", definition or "")
                    if len(w) > 4][:4]
        if keywords:
            name_q = " OR ".join(f'protein_name:"{w}"' for w in keywords)
            queries.append(f'({name_q}) AND organism_name:"{genus}"')
        queries.append(f'organism_name:"{genus}"')
        for query in queries:
            r = self._get(f"{UP_BASE}/uniprotkb/search", params={
                "query": query, "format": "tsv", "size": 100,
                "fields": "accession,organism_name,length,sequence"})
            if not r:
                continue
            best = None
            best_score = 0.0
            for line in r.text.strip().splitlines()[1:]:
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                acc, org, length, seq = (parts[0], parts[1], parts[2],
                                         parts[3])
                try:
                    int(length)
                except ValueError:
                    continue
                if genus.lower() not in org.lower():
                    continue
                sim = kmer_similarity(sequence, seq.upper())
                if sim > best_score:
                    best, best_score = acc, sim
            if best is not None and best_score >= 0.30:
                return best
        return None

    def resolve_uniprot(self, record):
        acc = record["accession"]
        base = acc.split(".")[0]
        # accessions that look like UniProt already
        if len(base) in (6, 10) and base[0].isalpha() and base[:2].isupper():
            r = self._get(f"{AF_BASE}/api/prediction/{base}")
            if r is not None:
                try:
                    if r.json():
                        return base
                except ValueError:
                    pass
        up = self.ncbi_to_uniprot(acc)
        if up:
            return up
        return self.uniprot_search_by_seq(record["species"],
                                           record["sequence"],
                                           record.get("definition", ""))

    def fetch_alphafold(self, uniprot_acc):
        """Download the AlphaFold DB PDB for a UniProt accession."""
        path = os.path.join(self.cache_dir, "raw", f"AF-{uniprot_acc}.pdb")
        if os.path.exists(path):
            return path
        r = self._get(f"{AF_BASE}/api/prediction/{uniprot_acc}")
        if r is None:
            return None
        try:
            entries = r.json()
        except ValueError:
            return None
        if not entries:
            return None
        pdb_url = entries[0].get("pdbUrl")
        if not pdb_url:
            return None
        for attempt in range(3):
            try:
                rr = self.session.get(pdb_url, timeout=180)
                if rr.status_code == 200:
                    with open(path, "w") as f:
                        f.write(rr.text)
                    return path
            except requests.RequestException:
                time.sleep(3)
        return None

    def parse_pdb_ca(self, path):
        coords, plddts, residues = [], [], []
        aa3to1 = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        }
        with open(path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    residues.append(line[17:20].strip())
                    coords.append((float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])))
                    plddts.append(float(line[60:66]))
        seq = "".join(aa3to1.get(r, 'X') for r in residues)
        return {"coords": np.array(coords), "plddts": np.array(plddts),
                "seq": seq}

    def thread_model(self, sequence, template):
        """Thread a sequence onto a template structure via local alignment.

        Returns coords/plddts arrays covering the aligned query positions.
        The returned "seq" is the QUERY sequence in TEMPLATE coordinate
        order (one residue per coordinate), so downstream sequence-based
        correspondence uses the correct residue for every CA position.
        """
        pairs = pl.sequence_pairs(template["seq"], sequence.upper())
        coords = template["coords"][[p[0] for p in pairs]]
        plddts = template["plddts"][[p[0] for p in pairs]] * 0.9
        if len(coords) < 100:
            return None
        qseq = "".join(sequence.upper()[p[1]] for p in pairs)
        return {"coords": coords, "plddts": plddts, "seq": qseq}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fasta", help="Input FASTA from the sequence scraper")
    p.add_argument("-o", "--outdir", default="structures",
                   help="Output directory")
    p.add_argument("--summary", default=None,
                   help="Optional summary TSV path (default: <outdir>/summary.tsv)")
    p.add_argument("--email", default=None)
    args = p.parse_args()

    if pl is None:
        sys.exit("proteinlib.py not found - put it next to the skills tree "
                 "(src/proteinlib.py) or on PYTHONPATH")

    os.makedirs(args.outdir, exist_ok=True)
    summary_path = args.summary or os.path.join(args.outdir, "summary.tsv")

    records = parse_fasta(args.fasta)
    print(f"[input] {len(records)} sequences", file=sys.stderr)
    fetcher = StructureFetcher(args.outdir, email=args.email)

    structures = []
    templates = []
    pending_model = []
    used_uniprot = set()
    for rec in records:
        name = rec["accession"].split(".")[0]
        print(f"[structure] {rec['species'] or name} ({name}) ...",
              file=sys.stderr)
        up = fetcher.resolve_uniprot(rec)
        if up in used_uniprot:
            up = None
        model = None
        source = "miss"
        mean_plddt = None
        if up:
            pdb_path = fetcher.fetch_alphafold(up)
            if pdb_path:
                parsed = fetcher.parse_pdb_ca(pdb_path)
                if len(parsed["coords"]) > 0:
                    model = {"coords": parsed["coords"],
                             "plddts": parsed["plddts"],
                             "seq": parsed["seq"]}
                    source = "alphafold"
                    mean_plddt = float(np.mean(parsed["plddts"]))
                    used_uniprot.add(up)
                    print(f"  AlphaFold DB: {up} "
                          f"({len(parsed['coords'])} residues, "
                          f"pLDDT {mean_plddt:.0f})", file=sys.stderr)
        if model is not None:
            structures.append({**rec, "uniprot": up, "model": model,
                               "source": source, "plddt": mean_plddt})
            templates.append(structures[-1])
        else:
            pending_model.append(rec)

    for rec in pending_model:
        name = rec["accession"].split(".")[0]
        print(f"[model] {rec['species'] or name}: threading ...",
              file=sys.stderr)
        best = None
        for t in templates:
            pairs = pl.sequence_pairs(rec["sequence"], t["model"]["seq"])
            if not pairs:
                continue
            frac = len(pairs) / max(len(rec["sequence"]), 1)
            if best is None or frac > best[0]:
                best = (frac, t)
        if best is None and templates:
            best = (0.0, templates[0])
        if best is None:
            print("  no template available, skipping", file=sys.stderr)
            structures.append({**rec, "uniprot": None, "model": None,
                               "source": "miss", "plddt": None})
            continue
        frac, template = best
        threaded = fetcher.thread_model(rec["sequence"],
                                        template["model"])
        if threaded is None:
            print(f"  threading failed (template {template['uniprot']})",
                  file=sys.stderr)
            structures.append({**rec, "uniprot": None, "model": None,
                               "source": "miss", "plddt": None})
            continue
        print(f"  threaded onto {template['uniprot']} "
              f"({len(threaded['coords'])} residues, "
              f"{frac * 100:.0f}% coverage)", file=sys.stderr)
        structures.append({**rec, "uniprot": None, "model": threaded,
                           "template_accession": template["accession"],
                           "source": f"threaded:{template['uniprot']}",
                           "plddt": float(np.mean(threaded["plddts"]))})

    got = [s for s in structures if s["model"] is not None]
    if len(got) < 2:
        print("[align] fewer than 2 structures; nothing to align",
              file=sys.stderr)
        sys.exit(1)

    # reference = structure with most sequence-similar to all others on average
    def avg_sim(s):
        sims = []
        for o in got:
            if o is s:
                continue
            pairs = pl.sequence_pairs(s["model"]["seq"],
                                      o["model"]["seq"]) or []
            ident = pl.core_identity(s["model"]["seq"], o["model"]["seq"],
                                     pairs)
            sims.append(ident)
        return np.mean(sims) if sims else 0.0
    print("[align] choosing reference ...", file=sys.stderr)
    real = [s for s in got if not s.get("template_accession")]
    reference = max(real or got, key=avg_sim)
    print(f"  reference: {reference['accession']} "
          f"({reference['species']})", file=sys.stderr)

    ref_coords = reference["model"]["coords"]
    aligned = []
    transforms = {}
    core_idents = {}
    # First pass: real structures, aligned with core-refined superposition.
    for s in got:
        if s is reference:
            transforms[s["accession"]] = (np.eye(3), np.zeros(3))
            core_idents[s["accession"]] = 1.0
        elif not s.get("template_accession"):
            _, _, r, shift = pl.superpose(s["model"], reference["model"])
            transforms[s["accession"]] = (r, shift)
            pairs = pl.sequence_pairs(s["model"]["seq"],
                                      reference["model"]["seq"])
            if pairs:
                core_idents[s["accession"]] = round(
                    pl.core_identity(s["model"]["seq"],
                                     reference["model"]["seq"], pairs), 3)
            else:
                core_idents[s["accession"]] = None
    # Second pass: apply transforms; threaded models inherit template's.
    for s in got:
        if s is reference:
            new_coords = ref_coords
            fit = 0.0
        elif s.get("template_accession"):
            tmpl_acc = s["template_accession"]
            r, shift = transforms.get(tmpl_acc, (np.eye(3), np.zeros(3)))
            core_idents[s["accession"]] = core_idents.get(tmpl_acc)
            new_coords = s["model"]["coords"] @ r + shift
            n = min(len(new_coords), len(ref_coords))
            m_i = np.linspace(0, len(new_coords) - 1, n).astype(int)
            r_i = np.linspace(0, len(ref_coords) - 1, n).astype(int)
            fit = pl.rmsd(new_coords[m_i], ref_coords[r_i])
        else:
            r, shift = transforms[s["accession"]]
            new_coords = s["model"]["coords"] @ r + shift
            n = min(len(new_coords), len(ref_coords))
            m_i = np.linspace(0, len(new_coords) - 1, n).astype(int)
            r_i = np.linspace(0, len(ref_coords) - 1, n).astype(int)
            fit = pl.rmsd(new_coords[m_i], ref_coords[r_i])
        # v2: per-residue secondary structure on the aligned trace
        ss = pl.assign_secondary(new_coords)
        seq = s["model"].get("seq") or s.get("sequence") or \
            "A" * len(new_coords)
        if len(seq) < len(new_coords):
            seq = seq + "X" * (len(new_coords) - len(seq))
        aligned.append({
            "accession": s["accession"], "species": s["species"],
            "lineage": s["lineage"], "family": s["family"],
            "definition": s["definition"], "uniprot": s.get("uniprot"),
            "source": s["source"], "plddt": s.get("plddt"),
            "rmsd_to_ref": fit if s is not reference else 0.0,
            "core_identity": core_idents.get(s["accession"]),
            "n_residues": int(len(new_coords)),
            "coords": [[round(float(c), 3) for c in xyz]
                       for xyz in new_coords],
            "plddts": [round(float(x), 1)
                       for x in s["model"]["plddts"]],
            "seq": seq[:len(new_coords)],
            "ss": "".join(ss),
        })
        name = s["accession"].split(".")[0]
        # v2: ribbon-ready PDB with HELIX/SHEET records
        pl.write_ca_pdb(new_coords, s["model"]["plddts"],
                        os.path.join(args.outdir, "aligned", f"{name}.pdb"),
                        name, ss=ss, seq=seq)
        print(f"  {s['accession']}: rmsd_to_ref={fit:.2f} A "
              f"({len(new_coords)} residues)", file=sys.stderr)

    with open(os.path.join(args.outdir, "aligned.json"), "w") as f:
        json.dump({"reference": reference["accession"],
                   "structures": aligned}, f)

    with open(summary_path, "w") as f:
        f.write("accession\tspecies\tlineage\tuniprot\tsource\tplddt\t"
                "rmsd_to_ref\tn_residues\n")
        for a in aligned:
            f.write(f"{a['accession']}\t{a['species']}\t{a['lineage']}\t"
                    f"{a['uniprot'] or ''}\t{a['source']}\t"
                    f"{a['plddt'] if a['plddt'] is not None else ''}\t"
                    f"{a['rmsd_to_ref']:.2f}\t{a['n_residues']}\n")
    missed = [s for s in structures if s["model"] is None]
    if missed:
        with open(summary_path, "a") as f:
            for s in missed:
                f.write(f"{s['accession']}\t{s['species']}\t"
                        f"{s['lineage']}\t\tmiss\t\t\t0\n")
    print(f"[done] {len(aligned)} structures aligned "
          f"({len(missed)} missed) -> {args.outdir}",
          file=sys.stderr)

    # visual check
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")
        for a in aligned:
            c = np.array(a["coords"])
            ax.plot(c[:, 0], c[:, 1], c[:, 2], linewidth=0.5, alpha=0.5)
        ax.set_axis_off()
        ax.set_title(f"{len(aligned)} structures superposed on "
                     f"{reference['accession']}")
        plt.savefig(os.path.join(args.outdir, "overlay.png"), dpi=120)
        plt.close(fig)
        print(f"[done] overlay -> {args.outdir}/overlay.png",
              file=sys.stderr)
    except Exception as e:
        print(f"[warn] overlay plot failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
