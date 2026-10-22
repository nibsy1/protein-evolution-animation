import os
import sys
import time
import requests

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, ".."),
              os.path.join(_HERE, "..", "..", "src"),
              os.path.join(_HERE, "..", "..", "..", "src")):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, os.path.abspath(_cand))
try:
    import proteinlib as pl
except ImportError:
    pl = None

from config import Config

AF_BASE = "https://alphafold.ebi.ac.uk"
UP_BASE = "https://rest.uniprot.org"
OUT_DIR = os.path.join("output", "structures")


class StructureAgent:
    """Fetches AlphaFold DB structures for homologs; builds homology models
    for species without a published AlphaFold entry."""

    def __init__(self, state):
        self.state = state
        os.makedirs(OUT_DIR, exist_ok=True)
        self.session = requests.Session()

    def _get(self, url, params=None, data=None, headers=None, timeout=60):
        for attempt in range(3):
            try:
                r = self.session.get(url, params=params, headers=headers, timeout=timeout)
                if r.status_code == 200:
                    return r
            except requests.RequestException:
                pass
            time.sleep(2 * (attempt + 1))
        return None

    def map_to_uniprot(self, accessions):
        """Map NCBI protein accessions to UniProt accessions."""
        r = requests.post(f"{UP_BASE}/idmapping/run", timeout=30,
                          data={"from": "RefSeq_Protein", "to": "UniProtKB",
                                "ids": ",".join(accessions)})
        if r.status_code != 200:
            return {}
        job_id = r.json()["jobId"]
        for _ in range(30):
            time.sleep(3)
            res = self._get(
                f"{UP_BASE}/idmapping/uniprotkb/results/stream/{job_id}",
                params={"format": "tsv", "fields": "accession,organism_name"},
            )
            if res and "From" in res.text:
                mapping = {}
                for line in res.text.strip().splitlines()[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        mapping[parts[0]] = parts[1]
                return mapping
        return {}

    def search_uniprot_fallback(self, organism, sequence, exclude=()):
        """Find a UniProt entry for an organism by protein name and length."""
        genus = organism.split()[0]
        query = (f'protein_name:"Calcium-transporting ATPase" AND '
                 f'(organism_name:"{organism}" OR organism_name:"{genus}")')
        res = self._get(f"{UP_BASE}/uniprotkb/search", params={
            "query": query, "format": "tsv",
            "fields": "accession,organism_name,length", "size": 20,
        })
        if not res:
            return None
        candidates = []
        for line in res.text.strip().splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            acc, org, length = parts[0], parts[1], parts[2]
            try:
                diff = abs(int(length) - len(sequence))
            except ValueError:
                continue
            if genus.lower() not in org.lower():
                continue
            candidates.append((diff, acc))
        candidates.sort()
        for _, acc in candidates:
            if acc not in exclude:
                return acc
        return None

    def fetch_alphafold(self, uniprot_acc):
        res = self._get(f"{AF_BASE}/api/prediction/{uniprot_acc}")
        if not res:
            return None
        try:
            entries = res.json()
        except ValueError:
            return None
        return entries[0] if entries else None

    def download_pdb(self, uniprot_acc):
        path = os.path.join(OUT_DIR, f"AF-{uniprot_acc}.pdb")
        if os.path.exists(path):
            return path
        for version in ["model_v6", "model_v4", "model_v3", "model_v2"]:
            url = f"{AF_BASE}/files/AF-{uniprot_acc}-F1-{version}.pdb"
            try:
                r = self.session.get(url, timeout=120)
                if r.status_code == 200:
                    with open(path, "w") as f:
                        f.write(r.text)
                    return path
            except requests.RequestException:
                continue
        return None

    def parse_pdb_ca(self, path):
        coords, plddts, residues = [], [], []
        with open(path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    residues.append(line[17:20].strip())
                    coords.append((float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])))
                    plddts.append(float(line[60:66]))
        return {"coords": coords, "plddts": plddts, "residues": residues}

    def build_homology_model(self, sequence, template_pdb):
        """Thread a sequence onto the closest template structure (CA trace).

        v2: uses a LOCAL BLOSUM62 alignment (via proteinlib) - global
        alignments are corrupted by unrelated terminal extensions on
        distant homologs.
        """
        template = self.parse_pdb_ca(template_pdb)
        if pl is not None:
            pairs = pl.sequence_pairs(template["seq"], sequence.upper())
        else:
            from Bio import Align
            from Bio.Align import PairwiseAligner
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aligner.open_gap_score = -11
            aligner.extend_gap_score = -1
            try:
                alignment = aligner.align(template["seq"],
                                          sequence.upper())[0]
            except Exception:
                return None
            pairs = []
            t_idx, q_idx = 0, 0
            for t_range, q_range in zip(alignment.aligned[0],
                                        alignment.aligned[1]):
                while t_idx < t_range[0]:
                    t_idx += 1
                while q_idx < q_range[0]:
                    q_idx += 1
                while t_idx < t_range[1] and q_idx < q_range[1]:
                    pairs.append((t_idx, q_idx))
                    t_idx += 1
                    q_idx += 1
        t_coords = np.asarray(template["coords"], dtype=float)
        t_plddt = np.asarray(template["plddts"], dtype=float)
        t_idx = [p[0] for p in pairs if p[0] < len(t_coords)]
        coords = t_coords[t_idx]
        plddts = t_plddt[t_idx] * 0.9
        if len(coords) < 100:
            return None
        return {"coords": coords, "plddts": plddts, "mode": "homology"}

    def run(self):
        homologs = self.state.get_homologs()
        structures = self.state.data.get("structures", [])
        if structures:
            print(f"   [structures] {len(structures)} already saved, skipping")
            return structures

        accessions = [h["accession"] for h in homologs]
        print(f"   [structures] Mapping {len(accessions)} accessions to UniProt...")
        mapping = self.map_to_uniprot(accessions)
        print(f"   [structures] Direct mappings: {len(mapping)}")

        results = []
        template_cache = {}
        used_uniprot = set()
        needs_model = []
        for homolog in homologs:
            lineage = homolog["lineage"]
            up_acc = mapping.get(homolog["accession"])
            if not up_acc and homolog["accession"].split(".")[0].startswith(("Q", "P", "O")):
                up_acc = homolog["accession"].split(".")[0]
            if not up_acc:
                up_acc = self.search_uniprot_fallback(homolog["organism"],
                                                      homolog["sequence"],
                                                      exclude=used_uniprot)
            if up_acc in used_uniprot:
                up_acc = None
            entry = self.fetch_alphafold(up_acc) if up_acc else None
            if entry:
                pdb_path = self.download_pdb(up_acc)
                if pdb_path:
                    parsed = self.parse_pdb_ca(pdb_path)
                    parsed["mode"] = "alphafold_db"
                    mean_plddt = entry.get("globalMetricValue")
                    print(f"   \u2714 {lineage}: {up_acc} ({homolog['organism']}) "
                          f"AlphaFold DB, pLDDT {mean_plddt:.1f}")
                    results.append({
                        "lineage": lineage,
                        "uniprot": up_acc,
                        "organism": homolog["organism"],
                        "ncbi_accession": homolog["accession"],
                        "mean_plddt": mean_plddt,
                        "mode": "alphafold_db",
                        "coords": parsed["coords"],
                        "plddts": parsed["plddts"],
                        "seq": parsed.get("seq", homolog["sequence"]),
                    })
                    used_uniprot.add(up_acc)
                    template_cache[lineage] = results[-1]
                    continue
            needs_model.append((homolog, up_acc))
            if up_acc:
                print(f"   \u26a0 {lineage}: {up_acc} has no AlphaFold DB entry, "
                      f"queueing homology model")
            else:
                print(f"   \u26a0 {lineage}: no distinct UniProt entry for "
                      f"{homolog['accession']} ({homolog['organism']}), "
                      f"queueing homology model")

        for homolog, up_acc in needs_model:
            lineage = homolog["lineage"]
            template = self._best_template(lineage, template_cache, results)
            if not template:
                print(f"   \u2718 {lineage}: no template available")
                continue
            template_pdb = (os.path.join(OUT_DIR, f"AF-{template['uniprot']}.pdb")
                            if template.get("mode") == "alphafold_db"
                            else None)
            if not template_pdb or not os.path.exists(template_pdb):
                print(f"   \u2718 {lineage}: template PDB unavailable")
                continue
            model = self.build_homology_model(homolog["sequence"], template_pdb)
            if not model or len(model["coords"]) < 100:
                print(f"   \u2718 {lineage}: homology model failed")
                continue
            print(f"   \u2714 {lineage}: homology model for {homolog['organism']} "
                  f"({len(model['coords'])} residues threaded onto "
                  f"{template['uniprot']})")
            results.append({
                "lineage": lineage,
                "uniprot": up_acc or homolog["accession"],
                "organism": homolog["organism"],
                "ncbi_accession": homolog["accession"],
                "mean_plddt": sum(model["plddts"]) / len(model["plddts"]),
                "mode": "homology",
                "coords": model["coords"],
                "plddts": model["plddts"],
                "seq": homolog["sequence"],
            })
            template_cache[lineage] = results[-1]

        self.state.add_structures(results)
        self.state.update_step(3)
        return results

    def _best_template(self, lineage, template_cache, results):
        order = Config.EVOLUTIONARY_STAGES
        idx = order.index(lineage) if lineage in order else 0
        for offset in range(1, len(order)):
            for candidate_idx in (idx - offset, idx + offset):
                if 0 <= candidate_idx < len(order):
                    key = order[candidate_idx]
                    if key in template_cache:
                        return template_cache[key]
        for res in results:
            return res
        return None
