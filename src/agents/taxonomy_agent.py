import json
import requests

from config import Config


class TaxonomyAgent:
    """Plans organism coverage across evolutionary lineages.

    Tries Mistral to select species; falls back to a curated list per lineage.
    """

    FALLBACK = {
        "Chlorophyta": [
            {"name": "Chlamydomonas reinhardtii", "taxon": 3055},
            {"name": "Ostreococcus tauri", "taxon": 70416},
        ],
        "Bryophyta": [
            {"name": "Physcomitrium patens", "taxon": 3218},
            {"name": "Marchantia polymorpha", "taxon": 3197},
        ],
        "Pteridophyta": [
            {"name": "Selaginella moellendorffii", "taxon": 158233},
            {"name": "Ceratopteris richardii", "taxon": 59132},
        ],
        "Pinophyta": [
            {"name": "Picea abies", "taxon": 69839},
            {"name": "Pinus taeda", "taxon": 3352},
        ],
        "Magnoliophyta": [
            {"name": "Arabidopsis thaliana", "taxon": 3702},
            {"name": "Oryza sativa", "taxon": 4530},
        ],
    }

    def __init__(self):
        self.api_key = None
        try:
            from dotenv import load_dotenv
            import os
            load_dotenv()
            self.api_key = os.getenv("MISTRAL_API_KEY")
        except ImportError:
            pass

    def plan_organisms(self):
        organisms = self._plan_with_mistral()
        if not organisms:
            organisms = self._plan_fallback()
        return organisms

    def _plan_with_mistral(self):
        if not self.api_key:
            return []
        headers = {"Authorization": f"Bearer {self.api_key}",
                   "Content-Type": "application/json"}
        prompt = (
            f"Select {Config.SPECIES_PER_LINEAGE} species per lineage for studying "
            f"the protein {Config.TARGET_PROTEIN}. Lineages: {Config.EVOLUTIONARY_STAGES}. "
            "Return ONLY a JSON array of objects with keys: name, lineage. "
            "Prefer species with sequenced genomes."
        )
        body = {
            "model": Config.MISTRAL_MODEL,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.2,
        }
        try:
            r = requests.post("https://api.mistral.ai/v1/chat/completions",
                              headers=headers, json=body, timeout=30)
            if r.status_code != 200:
                print(f"   [taxonomy] Mistral API {r.status_code}, using fallback list")
                return []
            text = r.json()["choices"][0]["message"]["content"]
            start = text.find("[")
            end = text.rfind("]")
            organisms = json.loads(text[start:end + 1])
            for org in organisms:
                org["lineage"] = org.get("lineage", "unknown")
            return organisms
        except Exception as e:
            print(f"   [taxonomy] Mistral failed ({e}), using fallback list")
            return []

    def _plan_fallback(self):
        organisms = []
        for lineage in Config.EVOLUTIONARY_STAGES:
            for org in self.FALLBACK.get(lineage, []):
                organisms.append({"name": org["name"], "lineage": lineage})
        return organisms
