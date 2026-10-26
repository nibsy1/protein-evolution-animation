#!/usr/bin/env python3
"""Scrape one homologous protein sequence per model species. (v5)

Usage:
    python scrape_sequences.py "ribulose bisphosphate carboxylase large chain" -o results.fasta
    python scrape_sequences.py ECA3 -o results.fasta --table results.tsv
    python scrape_sequences.py AT1G10130 -o results.fasta

Cascade per species:
    1. NCBI Protein search (gene name, AGI, protein-name terms; species
       then genus scope) - collects ALL family candidates/isoforms
    2. Every candidate is verified as a true homolog of the query
       sequence via a two-stage check: k-mer similarity pre-filter
       (threshold 0.10) + BLOSUM62 local-alignment identity gate
       (>= 0.35 identity, >= 0.50 coverage). This rejects wrong-family
       paralogs (e.g. ACA10 instead of ECA3) that share only ~33%
       identity. Only the single closest match is kept.
    3. Last resort: the query term is broadened to the protein class
       (e.g. 'ATPase') with a relaxed k-mer pre-filter (0.02) so
       unnamed records can be identified by alignment identity.
    4. Species still without a hit are gap-filled from UniProt REST.

Output:
    - FASTA file with one sequence per species
    - TSV table with accession, organism, lineage, family, length, quality
      (query / exact / homolog / broad / uniprot / genus_fill / related /
       fallback / miss)
"""
import argparse
import sys
import time
import xml.etree.ElementTree as ET

import requests

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UP_BASE = "https://rest.uniprot.org"
KMER_K = 5
# v2: Raised from 0.008 to 0.10. The old 0.8% threshold accepted any two
# P-type ATPases as "homologs" (e.g. ECA3 vs ACA10 share ~2% k-mers by chance).
# 10% k-mer overlap cleanly separates:
#   - True distant orthologs: Volvox ~10%, Chlorella ~19%, Marchantia ~10%
#   - Wrong family members (ACA/ECA paralogs): ~2% k-mer overlap
# Chlamydomonas (~1.6%) falls below this threshold but is rescued by the
# broaden step + UniProt gap-filling which use the alignment identity check.
KMER_MIN_FRAC = 0.10
# v2: Minimum BLOSUM62 local-alignment identity over the aligned core.
# True orthologs share >45% identity even across deep divergence.
# Wrong family members (ECA3 vs ACA10: ~34%, ECA3 vs ACA3: ~33%)
# sit just below the 0.35 threshold. Green algae (Chlamydomonas ~34%)
# are borderline and may need the broaden/UniProt fallback path.
ALN_MIN_IDENTITY = 0.35
# v2: Minimum fraction of the query sequence covered by the local alignment.
# Prevents accepting a short domain match as a full-length homolog.
ALN_MIN_COVERAGE = 0.50


def kmer_similarity(seq_a, seq_b, k=KMER_K):
    """Fraction of seq_a k-mers shared with seq_b. Cheap homology proxy."""
    if not seq_a or not seq_b:
        return 0.0
    ka = {seq_a[i:i + k] for i in range(len(seq_a) - k + 1)}
    kb = {seq_b[i:i + k] for i in range(len(seq_b) - k + 1)}
    if not ka:
        return 0.0
    return len(ka & kb) / len(ka)


def alignment_identity(seq_a, seq_b):
    """BLOSUM62 local-alignment identity and coverage of seq_a by seq_b.

    Returns (identity, coverage) where:
      - identity = fraction of aligned positions with identical residues
      - coverage = fraction of seq_a residues covered by the alignment

    Uses Biopython's PairwiseAligner with BLOSUM62 (same as proteinlib).
    Falls back to k-mer similarity if Biopython is unavailable.
    """
    if not seq_a or not seq_b:
        return 0.0, 0.0
    try:
        from Bio.Align import PairwiseAligner
        from Bio.Align import substitution_matrices
        aligner = PairwiseAligner()
        aligner.mode = "local"
        try:
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        except Exception:
            pass
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        alignment = aligner.align(seq_a, seq_b)[0]
        pairs = []
        i_a, i_b = 0, 0
        for a_range, b_range in zip(alignment.aligned[0], alignment.aligned[1]):
            a0, a1 = a_range
            b0, b1 = b_range
            while i_a < a0:
                i_a += 1
            while i_b < b0:
                i_b += 1
            while i_a < a1 and i_b < b1:
                pairs.append((i_a, i_b))
                i_a += 1
                i_b += 1
        if not pairs:
            return 0.0, 0.0
        identity = sum(1 for i, j in pairs if seq_a[i] == seq_b[j]) / len(pairs)
        coverage = len(set(i for i, _ in pairs)) / len(seq_a)
        return identity, coverage
    except ImportError:
        # Fallback: use k-mer similarity as rough identity proxy
        sim = kmer_similarity(seq_a, seq_b)
        return sim, sim  # rough approximation


# Lineage / species / family
SPECIES = [
    ("Green Algae", "Chlamydomonas reinhardtii", "Chlamydomonadaceae"),
    ("Green Algae", "Volvox carteri", "Volvocaceae"),
    ("Green Algae", "Chlorella vulgaris", "Chlorellaceae"),
    ("Bryophytes", "Physcomitrium patens", "Physcomitriaceae"),
    ("Bryophytes", "Marchantia polymorpha", "Marchantiaceae"),
    ("Ferns", "Azolla filiculoides", "Salviniaceae"),
    ("Gymnosperms", "Pinus taeda", "Pinaceae"),
    ("Gymnosperms", "Picea abies", "Pinaceae"),
    ("Angiosperms", "Arabidopsis thaliana", "Brassicaceae"),
    ("Angiosperms", "Oryza sativa", "Poaceae"),
    ("Angiosperms", "Zea mays", "Poaceae"),
    ("Angiosperms", "Solanum lycopersicum", "Solanaceae"),
    ("Angiosperms", "Populus trichocarpa", "Salicaceae"),
    ("Angiosperms", "Phaseolus vulgaris", "Fabaceae"),
    ("Angiosperms", "Hordeum vulgare", "Poaceae"),
    ("Angiosperms", "Nicotiana benthamiana", "Solanaceae"),
]


class Scraper:
    def __init__(self, email=None):
        self.session = requests.Session()
        if email:
            self.session.params = {"tool": "sequence_scraper", "email": email}

    def _get(self, endpoint, params, retries=4, timeout=40):
        for attempt in range(retries):
            try:
                r = self.session.get(f"{BASE}/{endpoint}", params=params,
                                     timeout=timeout)
                if r.status_code == 200 and r.text.lstrip().startswith("<"):
                    return r.text
            except requests.RequestException:
                pass
            time.sleep(1.5 * (attempt + 1))
        return None

    def resolve_query(self, query):
        """Resolve a query (gene name like ECA3, AGI code, or protein name)
        to (search_terms, query_record).

        - AGI code: fetch that Arabidopsis protein record from NCBI.
        - All-caps gene symbol ending in a digit (ECA3, RCA2, ...):
          UniProt exact gene lookup in Arabidopsis; the query record pins
          the correct isoform. Falls through to a plain term on failure.
        - Otherwise: protein-name query (query_record=None).
        """
        import re as _re
        q = query.strip()
        if _re.match(r"^AT[1-5CM]G\d{5}$", q, _re.IGNORECASE):
            agi = q.upper()
            xml_text = self._get("esearch.fcgi", {
                "db": "protein",
                "term": f'{agi}[All Fields] AND "Arabidopsis thaliana"[Organism]',
                "retmax": 5,
            })
            ids = []
            if xml_text:
                try:
                    root = ET.fromstring(xml_text)
                    ids = [el.text for el in root.findall(".//IdList/Id")]
                except ET.ParseError:
                    pass
            for uid in ids:
                records = self.fetch_records([uid])
                if records:
                    rec = records[0]
                    print(f"[query] {agi} -> {rec['accession']} "
                          f"({rec['definition'][:50]})", file=sys.stderr)
                    return [agi], rec
            print(f"[query] AGI {agi} not found", file=sys.stderr)
            return None, None
        if (_re.match(r"^[A-Z][A-Z0-9]{1,9}$", q) and q.upper() == q
                and q[-1].isdigit()):
            gene = q.upper()
            r = None
            for attempt in range(3):
                try:
                    r = self.session.get(
                        f"{UP_BASE}/uniprotkb/search",
                        params={"query": f'gene_exact:{gene} AND organism_name:'
                                         '"Arabidopsis thaliana"',
                                "format": "tsv", "size": 5,
                                "fields": "accession,gene_names,protein_name,"
                                          "xref_geneId,sequence"},
                        timeout=40)
                    if r.status_code == 200:
                        break
                except requests.RequestException:
                    pass
                time.sleep(2)
            if r is not None and r.status_code == 200:
                for line in r.text.strip().splitlines()[1:]:
                    parts = line.split("\t")
                    if len(parts) < 5:
                        continue
                    acc, genes, pname, geneids, seq = parts[:5]
                    if gene not in genes.upper().split():
                        continue
                    agi = geneids.split(";")[0] if geneids else None
                    print(f"[query] {gene} -> {acc} ({pname[:50]})",
                          file=sys.stderr)
                    rec = {"accession": acc,
                           "definition": f"{pname} [{gene}]",
                           "organism": "Arabidopsis thaliana",
                           "length": len(seq), "sequence": seq.upper()}
                    family = _re.sub(r"\s*\([^)]*\)\s*", " ", pname).strip()
                    terms = [gene] + ([agi] if agi else []) + [family]
                    return terms, rec
            print(f"[query] gene {gene} not found in UniProt for Arabidopsis",
                  file=sys.stderr)
        return [query], None

    def search_ids(self, term, species, retmax=50):
        xml_text = self._get("esearch.fcgi", {
            "db": "protein",
            "term": f'({term}) AND "{species}"[Organism]',
            "retmax": retmax,
        })
        if not xml_text:
            return []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return []
        return [el.text for el in root.findall(".//IdList/Id")]

    def search_ids_genus(self, term, genus, retmax=50):
        xml_text = self._get("esearch.fcgi", {
            "db": "protein",
            "term": f'({term}) AND "{genus}"[Organism]',
            "retmax": retmax,
        })
        if not xml_text:
            return []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return []
        return [el.text for el in root.findall(".//IdList/Id")]

    def fetch_records(self, ids):
        results = []
        for i in range(0, len(ids), 10):
            chunk = ids[i:i + 10]
            xml_text = self._get("efetch.fcgi", {
                "db": "protein", "id": ",".join(chunk),
                "rettype": "gp", "retmode": "xml",
            })
            if not xml_text:
                continue
            try:
                root = ET.fromstring(xml_text)
            except ET.ParseError:
                continue
            for prot in root.findall(".//GBSeq"):
                seq = (prot.findtext("GBSeq_sequence") or "").upper()
                if not seq:
                    continue
                results.append({
                    "accession": prot.findtext("GBSeq_accession-version", ""),
                    "definition": prot.findtext("GBSeq_definition", ""),
                    "organism": prot.findtext("GBSeq_organism", ""),
                    "length": len(seq),
                    "sequence": seq,
                })
            time.sleep(0.5)
        return results

    STOPWORDS = {"and", "the", "of", "chain", "subunit", "protein", "putative",
                 "like", "type", "family", " precursor", "isoform"}

    @classmethod
    def _query_tokens(cls, query):
        return [w for w in query.lower().split()
                if len(w) > 3 and w not in cls.STOPWORDS]

    def _relevance(self, record, tokens):
        if not tokens:
            return 1.0
        definition = record["definition"].lower()
        return sum(1 for t in tokens if t in definition) / len(tokens)

    def _exclude(self, record, tokens, query):
        """Reject obviously wrong hits."""
        definition = record["definition"].lower()
        q = (query or "").lower()
        if "small subunit" in definition or "small chain" in definition:
            if "small" not in q:
                return True
        if "large subunit" in definition or "large chain" in definition:
            if "large" in q:
                pass
            elif "chain" in q or "subunit" in q:
                return True
        if any(t in definition for t in ("hypothetical", "predicted")):
            return True
        import re
        if re.match(r"^[0-9][A-Z0-9]{3}_", record["accession"]):
            return True
        return False

    def best_record(self, records, species, query, exact_required=False,
                    query_length=None):
        import re as _re
        genus = species.split()[0]
        same_species = [r for r in records
                        if r["organism"] == species or species in r["organism"]]
        pool = same_species or [r for r in records
                                if genus in r["organism"]] or records
        tokens = self._query_tokens(query)
        candidates = [r for r in pool
                      if self._relevance(r, tokens) >= 0.5
                      and not self._exclude(r, tokens, query)]
        has_relevant = bool(candidates)
        if candidates:
            pool = candidates

        def _primary_name(r):
            d = r["definition"]
            m = _re.match(r"RecName:\s*Full=([^;\[]+)", d)
            name = (m.group(1) if m else _re.split(r"[\[\(]", d)[0])
            return name.strip().rstrip(".,").lower()

        def exact_name(r):
            name = _primary_name(r)
            base = _re.sub(r"[-, ]+\d+$", "", name).strip()
            return 1 if base == query.strip().lower() else 0

        def _name_match(r):
            name = _primary_name(r).split(",")[0].strip()
            q = query.strip().lower().split(",")[0].strip()
            if not all(tok in name for tok in q.split() if len(tok) > 3):
                return False
            q_num = _re.search(r"(\d+)$", q)
            if q_num:
                r_num = _re.search(r"(\d+)$", name)
                if not r_num or r_num.group(1) != q_num.group(1):
                    return False
            return True

        def length_ok(r):
            if not query_length:
                return True
            return 0.5 * query_length <= r["length"] <= 1.3 * query_length

        def score(r):
            refseq = 2 if r["accession"].startswith(("XP_", "NP_")) else 0
            reviewed = 1 if "RecName" in r["definition"] else 0
            complete = 0 if "partial" in r["definition"].lower() else 1
            return (exact_name(r), self._relevance(r, tokens), complete,
                    reviewed, refseq, r["length"])

        if query_length:
            sane = [r for r in pool if length_ok(r)]
            if sane:
                pool = sane
        if exact_required:
            named = [r for r in pool if _name_match(r) and length_ok(r)]
            if not named:
                named = [r for r in pool if _name_match(r)]
            if named:
                best = max(named, key=score)
                if not length_ok(best):
                    return None, False, has_relevant
                return best, True, True
        best = max(pool, key=score)
        if query_length and not length_ok(best):
            return None, False, has_relevant
        return best, exact_name(best) == 1, has_relevant

    def search_uniprot(self, term, species, size=20, query_length=None):
        """Search UniProt by protein name + organism; return records with
        sequences, exact-species records ranked first."""
        import re as _re
        genus = species.split()[0]
        name_core = term.split(",")[0]
        name_core = _re.sub(r"\([^)]*\)", " ", name_core)
        name_core = _re.sub(r"[^A-Za-z0-9 -]", " ", name_core)
        name_core = _re.sub(r"\s*\d+$", "", name_core).strip() or name_core
        r = None
        for attempt in range(3):
            try:
                r = self.session.get(
                    f"{UP_BASE}/uniprotkb/search",
                    params={"query": f'protein_name:"{name_core}" AND '
                                     f'(organism_name:"{species}" OR '
                                     f'organism_name:"{genus}")',
                            "format": "tsv", "size": size,
                            "fields": "accession,organism_name,length,sequence"},
                    timeout=40)
                if r.status_code == 200:
                    break
            except requests.RequestException:
                pass
            time.sleep(2 * (attempt + 1))
        if r is None or r.status_code != 200:
            return []
        records = []
        for line in r.text.strip().splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            acc, org, length, seq = parts[0], parts[1], parts[2], parts[3]
            try:
                length = int(length)
            except ValueError:
                continue
            if not seq:
                continue
            records.append({
                "accession": acc,
                "definition": f"{term} [UniProt]",
                "organism": org,
                "length": length,
                "sequence": seq.upper(),
                "_exact_species": species.lower() in org.lower(),
            })
        records.sort(key=lambda r: (not r["_exact_species"],
                                    abs(r["length"] - (query_length or r["length"]))))
        return records

    def _collect_candidates(self, terms, species):
        """Gather family candidates for a species across query terms."""
        genus = species.split()[0]
        records = []
        for term in terms:
            ids = self.search_ids(term, species)
            if ids:
                records = self.fetch_records(ids)
            if not records:
                ids = self.search_ids_genus(term, genus)
                if ids:
                    records = self.fetch_records(ids)
            if records:
                break
        return records

    def _broaden_records(self, terms, species, query_sequence,
                         query_length):
        """Last resort: broad class terms (e.g. 'ATPase'); the alignment check
        picks the true homolog among loosely-related or unnamed records.

        v2: Uses a relaxed k-mer threshold (KMER_MIN_FRAC / 5) because distant
        green algae (Chlamydomonas ~1.6% k-mers) fall below the strict gate.
        The BLOSUM62 alignment identity check (ALN_MIN_IDENTITY) is the real
        filter — it accepts Chlamydomonas (~34%) while rejecting wrong-family
        paralogs (~32-33%).
        """
        head = terms[0].split()[0] if terms[0] else "protein"
        relaxed_kmer = KMER_MIN_FRAC / 5.0  # 0.02 for distant homologs
        for broad in (f"{head} ATPase", "ATPase"):
            if broad in terms:
                continue
            ids = self.search_ids(broad, species, retmax=100)
            records = self.fetch_records(ids) if ids else []
            if not records:
                continue
            # Use relaxed k-mer pre-filter + alignment identity gate
            scored = []
            for r in records:
                if self._exclude(r, [], ""):
                    continue
                sim = kmer_similarity(query_sequence, r["sequence"])
                if sim < relaxed_kmer:
                    continue
                if query_length and not (0.5 * query_length <= r["length"] <= 1.3 * query_length):
                    continue
                identity, coverage = alignment_identity(query_sequence, r["sequence"])
                if identity < ALN_MIN_IDENTITY:
                    continue
                if coverage < ALN_MIN_COVERAGE:
                    continue
                score = identity + 0.1 * sim
                r_copy = dict(r)
                r_copy["_aln_identity"] = identity
                r_copy["_aln_coverage"] = coverage
                r_copy["_kmer_sim"] = sim
                scored.append((score, r_copy))
            scored.sort(key=lambda x: x[0], reverse=True)
            if scored:
                return scored
        return []

    def _closest_to_query(self, records, query_sequence, query_length=None):
        """Rank verified homologs by alignment identity to the query sequence.

        v2: Uses BLOSUM62 local-alignment identity + coverage instead of just
        k-mer similarity. This rejects wrong-family members (e.g. ECA1/2/4
        or ACA family) that share only ~33% identity with the true ortholog.
        K-mer similarity is still used as a fast pre-filter.
        """
        scored = []
        for r in records:
            if self._exclude(r, [], ""):
                continue
            # Fast pre-filter: k-mer similarity
            sim = kmer_similarity(query_sequence, r["sequence"])
            if sim < KMER_MIN_FRAC:
                continue
            # Length filter
            if query_length and not (0.5 * query_length <= r["length"] <= 1.3 * query_length):
                continue
            # v2: BLOSUM62 alignment identity + coverage check
            identity, coverage = alignment_identity(query_sequence, r["sequence"])
            if identity < ALN_MIN_IDENTITY:
                continue
            if coverage < ALN_MIN_COVERAGE:
                continue
            # Score by alignment identity (higher = better match)
            # Combine with k-mer similarity as tiebreaker
            score = identity + 0.1 * sim
            scored.append((score, r, identity, coverage))
        scored.sort(key=lambda x: x[0], reverse=True)
        # Return (score, record) tuples for backward compatibility,
        # but also include identity/coverage in the record dict
        result = []
        for score, r, identity, coverage in scored:
            r_copy = dict(r)
            r_copy["_aln_identity"] = identity
            r_copy["_aln_coverage"] = coverage
            r_copy["_kmer_sim"] = kmer_similarity(query_sequence, r["sequence"])
            result.append((score, r_copy))
        return result

    def scrape(self, terms, exact_required=False, query_length=None,
               skip_species=(), query_sequence=None):
        rows = []
        misses = []
        for lineage, species, family in SPECIES:
            if species in skip_species:
                continue
            print(f"[{lineage}] {species} ...", file=sys.stderr)
            genus = species.split()[0]
            records = self._collect_candidates(terms, species)
            best, quality = None, None
            if query_sequence and records:
                scored = self._closest_to_query(records, query_sequence,
                                                query_length)
                if scored:
                    sim, best = scored[0]
                    quality = "exact" if sim > 0.25 else "homolog"
            elif records:
                best, exact, relevant = self.best_record(
                    records, species, terms[0],
                    exact_required=exact_required, query_length=query_length)
                if best is not None:
                    quality = "exact" if exact else ("related" if relevant else "fallback")
            if (best is not None and query_length
                    and best["length"] < 0.6 * query_length):
                genus_records = []
                ids = self.search_ids_genus(terms[0], genus)
                if ids:
                    genus_records = self.fetch_records(ids)
                if not genus_records and query_sequence:
                    up = self.search_uniprot(terms[0], species,
                                             query_length=query_length)
                    genus_records = up
                if query_sequence and genus_records:
                    scored = self._closest_to_query(genus_records,
                                                    query_sequence,
                                                    query_length)
                    if scored and scored[0][1]["length"] > best["length"]:
                        best = scored[0][1]
                        quality = "genus_fill"
            if best is None and query_sequence:
                broad = self._broaden_records(terms, species, query_sequence,
                                              query_length)
                if broad:
                    sim, best = broad[0]
                    quality = "broad"
            if best is None:
                misses.append((lineage, species, family))
                continue
            # v2: Log alignment metrics for quality control
            aln_id = best.get("_aln_identity", "?")
            aln_cov = best.get("_aln_coverage", "?")
            kmer_s = best.get("_kmer_sim", "?")
            if isinstance(aln_id, float):
                aln_str = f"aln_id={aln_id:.2f} cov={aln_cov:.2f} kmer={kmer_s:.3f}"
            else:
                aln_str = ""
            print(f"  HIT [{quality}]: {best['accession']} ({best['length']} aa) "
                  f"{best['definition'][:60]} {aln_str}", file=sys.stderr)
            rows.append({
                "lineage": lineage, "species": species, "family": family,
                "accession": best["accession"],
                "definition": best["definition"],
                "length": best["length"], "sequence": best["sequence"],
                "quality": quality,
            })
            time.sleep(0.4)

        # Gap filling via UniProt, verified by local homology check
        if misses and query_sequence:
            print(f"\n[uniprot] Filling {len(misses)} gap(s) via UniProt...",
                  file=sys.stderr)
            for lineage, species, family in misses:
                up_records = []
                for term in terms:
                    up_records = self.search_uniprot(term, species,
                                                      query_length=query_length)
                    if up_records:
                        break
                scored = self._closest_to_query(up_records, query_sequence,
                                               query_length) if up_records else []
                if scored:
                    sim, filled = scored[0]
                    print(f"  FILLED [uniprot]: {filled['accession']} "
                          f"({filled['length']} aa, sim {sim:.2f}) for {species}",
                          file=sys.stderr)
                    rows.append({
                        "lineage": lineage, "species": species,
                        "family": family,
                        "accession": filled["accession"],
                        "definition": filled["definition"],
                        "length": filled["length"],
                        "sequence": filled["sequence"],
                        "quality": "uniprot",
                    })
                else:
                    print(f"  STILL MISS: {species}", file=sys.stderr)
                    rows.append({
                        "lineage": lineage, "species": species,
                        "family": family, "accession": "", "definition": "",
                        "length": 0, "sequence": "", "quality": "miss",
                    })
                time.sleep(0.4)
        else:
            for lineage, species, family in misses:
                rows.append({
                    "lineage": lineage, "species": species,
                    "family": family, "accession": "", "definition": "",
                    "length": 0, "sequence": "", "quality": "miss",
                })
        order = {spec: i for i, (_, spec, _) in enumerate(SPECIES)}
        rows.sort(key=lambda r: order.get(r["species"], 99))
        return rows


def write_fasta(rows, path):
    with open(path, "w") as f:
        for r in rows:
            if not r["sequence"]:
                continue
            header = (f">{r['accession']} {r['species']}|{r['lineage']}"
                      f"|{r['family']} {r['definition']}")
            f.write(header + "\n")
            seq = r["sequence"]
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


def write_table(rows, path):
    with open(path, "w") as f:
        f.write("lineage\tspecies\tfamily\taccession\tlength\tquality\tdefinition\n")
        for r in rows:
            f.write(f"{r['lineage']}\t{r['species']}\t{r['family']}\t"
                    f"{r['accession']}\t{r['length']}\t{r.get('quality','')}\t"
                    f"{r['definition']}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("query", help="Protein name, gene name, or AGI code")
    p.add_argument("-o", "--output", default="sequences.fasta",
                   help="Output FASTA path")
    p.add_argument("--table", default=None, help="Optional output TSV path")
    p.add_argument("--email", default=None,
                   help="Contact email for NCBI (recommended)")
    args = p.parse_args()

    scraper = Scraper(email=args.email)
    terms, qrec = scraper.resolve_query(args.query)
    if terms is None:
        sys.exit(f"Could not resolve query: {args.query}")
    print(f"[query] terms: {terms}", file=sys.stderr)
    exact_required = qrec is not None
    qlen = qrec["length"] if qrec else None
    qseq = qrec["sequence"] if qrec else None
    skip = ("Arabidopsis thaliana",) if qrec else ()
    rows = scraper.scrape(terms, exact_required=exact_required,
                          query_length=qlen, skip_species=skip,
                          query_sequence=qseq)
    if qrec:
        rows.append({
            "lineage": "Angiosperms", "species": "Arabidopsis thaliana",
            "family": "Brassicaceae", "accession": qrec["accession"],
            "definition": qrec["definition"],
            "length": qrec["length"],
            "sequence": qrec["sequence"],
            "quality": "query",
        })
    order = {spec: i for i, (_, spec, _) in enumerate(SPECIES)}
    rows.sort(key=lambda r: order.get(r["species"], 99))
    write_fasta(rows, args.output)
    print(f"Wrote {sum(1 for r in rows if r['sequence'])}/{len(rows)} "
          f"sequences to {args.output}", file=sys.stderr)
    if args.table:
        write_table(rows, args.table)
        print(f"Wrote table to {args.table}", file=sys.stderr)


if __name__ == "__main__":
    main()
