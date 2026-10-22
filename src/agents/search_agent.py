import time
import requests
import xml.etree.ElementTree as ET

from config import Config

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class SearchAgent:
    """Searches NCBI protein for homologs of the target protein."""

    def __init__(self, state):
        self.state = state
        self.session = requests.Session()

    def _get(self, endpoint, params, retries=4):
        for attempt in range(retries):
            try:
                r = self.session.get(f"{BASE}/{endpoint}", params=params, timeout=40)
                if r.status_code == 200 and r.text.lstrip().startswith("<"):
                    return r.text
            except requests.RequestException:
                pass
            time.sleep(1.5 * (attempt + 1))
        return None

    def _search(self, term, retmax):
        xml_text = self._get("esearch.fcgi", {
            "db": "protein", "term": term, "retmax": retmax, "sort": "Length",
        })
        if not xml_text:
            return []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return []
        return [id_el.text for id_el in root.findall(".//IdList/Id")]

    def _fetch(self, ids):
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
                accession = prot.findtext("GBSeq_accession-version", "")
                definition = prot.findtext("GBSeq_definition", "")
                organism = prot.findtext("GBSeq_organism", "")
                seq = (prot.findtext("GBSeq_sequence") or "").upper()
                if not seq:
                    continue
                results.append({
                    "accession": accession,
                    "definition": definition,
                    "organism": organism,
                    "sequence": seq,
                })
            time.sleep(0.5)
        return results

    def _species_queries(self, organism):
        target = Config.TARGET_PROTEIN
        name = organism["name"]
        genus = name.split()[0]
        scopes = [f'"{name}"', f'"{genus}"'] if genus != name else [f'"{name}"']
        queries = []
        for scope in scopes:
            queries.append(f'{target}[Gene Name] AND {scope}[Organism]')
            queries.append(f'calcium-transporting[All Fields] AND ATPase[All Fields] AND {scope}[Organism]')
            queries.append(f'ATPase[Title] AND calcium[Title] AND {scope}[Organism]')
        return queries

    def _lineage_queries(self, lineage):
        target = Config.TARGET_PROTEIN
        return [
            f'{target}[Gene Name] AND {lineage}[Organism]',
            f'calcium-transporting[All Fields] AND ATPase[All Fields] AND {lineage}[Organism]',
            f'{target}[Title] AND {lineage}[Organism]',
        ]

    def _match_found(self, record):
        definition = record["definition"].lower()
        return "atpase" in definition and ("calcium" in definition or "ca2" in definition
                                           or Config.TARGET_PROTEIN.lower() in definition)

    def _pick(self, records, exclude_species, wanted):
        picked = []
        for rec in sorted(records, key=lambda r: len(r["sequence"]), reverse=True):
            genus = rec["organism"].split()[0] if rec["organism"] else ""
            if genus in exclude_species:
                continue
            if len(rec["sequence"]) < 250:
                continue
            picked.append(rec)
            if len(picked) >= wanted:
                break
        return picked

    def execute_search(self, organisms):
        self.state.add_organisms(organisms)
        collected = []
        seen = set()
        taken_species = set()
        per_lineage = {}
        for organism in organisms:
            lineage = organism["lineage"]
            found = per_lineage.get(lineage, 0)
            if found >= Config.SPECIES_PER_LINEAGE:
                continue
            wanted = Config.SPECIES_PER_LINEAGE - found
            records = []
            for term in self._species_queries(organism):
                ids = self._search(term, retmax=10)
                if not ids:
                    time.sleep(0.4)
                    continue
                records = [r for r in self._fetch(ids) if self._match_found(r)]
                if records:
                    break
                time.sleep(0.4)
            source = organism["name"]
            if not records:
                for term in self._lineage_queries(lineage):
                    ids = self._search(term, retmax=20)
                    if not ids:
                        time.sleep(0.4)
                        continue
                    records = [r for r in self._fetch(ids) if self._match_found(r)]
                    if records:
                        source = f"{lineage} (lineage-wide search)"
                        break
                    time.sleep(0.4)
            picks = self._pick(records, taken_species, wanted)
            for rec in picks:
                rec["lineage"] = lineage
                rec["organism_planned"] = organism["name"]
                rec["source"] = source
                collected.append(rec)
                seen.add(rec["accession"])
                genus = rec["organism"].split()[0] if rec["organism"] else ""
                taken_species.add(genus)
                per_lineage[lineage] = per_lineage.get(lineage, 0) + 1
                print(f"   \u2714 {lineage}: {rec['accession']} ({rec['organism']}, "
                      f"{len(rec['sequence'])} aa)")
            if not picks:
                print(f"   \u2718 {lineage}: no homolog found for {organism['name']}")
            time.sleep(0.4)
        self.state.add_homologs(collected)
        return collected
