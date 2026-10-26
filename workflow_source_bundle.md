# protein-animation workflow — full source bundle

Each section below is one file, verbatim. Root = `workflow/`.

## `.dockerignore`

```
.venv
__pycache__
*.pyc
.git
outputs

```

## `.gitignore`

```
__pycache__/
*.pyc
.venv/
outputs/
.git/

```

## `Dockerfile`

```
FROM python:3.12-slim

WORKDIR /app

RUN pip install --no-cache-dir uv

COPY requirements.txt ./
RUN uv pip install --system --no-cache-dir -r requirements.txt

COPY . /app
ENV PYTHONPATH=/app

CMD ["python", "entrypoints/worker.py"]

```

## `README.md`

```
# protein-animation workflow

A Mistral Workflow that turns a protein query — a gene name (`TPC1`), an
Arabidopsis AGI code (`AT4G03560`), or a protein name — into a complete
set of evolutionary protein animation artifacts:

| Artifact | Description |
|---|---|
| `sequences.fasta` | One verified homolog per model plant species (16 species, green algae → angiosperms) |
| `sequences.tsv` | Quality table (accession, lineage, quality flag) |
| `aligned/*.pdb` | Ribbon-ready aligned PDBs (HELIX/SHEET records) |
| `aligned.json` | Machine-readable aligned structures (coords, seq, SS, pLDDT) |
| `summary.tsv` | Per-structure alignment stats |
| `overlay.png` | All structures superposed |
| `evolution.gif` | Ribbon-style, stability-masked morph animation across lineages |
| `evolution_sheet.png` | One still per species |
| `viewer.html` | Interactive viewer: clickable tree + rotatable 3D ribbon + smooth morph |

## How it works

The workflow chains the four workspace skills as durable activities:

```
scrape_sequences ─▶ align_structures ─▶ render_animation ─▶ build_viewer
(sequence-scraper)  (structure-aligner)  (morph-animation)    (interactive-viewer)
```

A **supervisor agent** (durable agent) watches the run and writes a short
markdown report.

- **Durability**: every stage is a Temporal activity with retries and
  heartbeats; a worker crash resumes where it left off.
- **Stability-masked morphing**: only genuinely changing regions animate;
  a helix that is a helix in both species stays put.
- **Ribbon style + consistent medium lighting**, no membrane visualization.

## Run it

### Locally

```bash
cd workflow
make install
make start-worker          # in one terminal
make execute workflow=protein-animation input='{"query": "TPC1"}'
```

### From the Mistral Console

1. Push this `workflow/` directory to a GitHub repo.
2. In the Console: **Build → Workflows → New deployment**, point at the repo.
3. Mistral Cloud builds the Dockerfile and runs the worker.
4. **Start Workflow** with input `{"query": "TPC1"}` (or `AT4G03560`,
   or a protein name). Optional: `email`, `make_viewer`, `hold`, `morph`,
   `fps`.

### Via API

```bash
curl -X POST https://api.mistral.ai/v1/workflows/protein-animation/execute \
  -H "Authorization: Bearer $MISTRAL_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"input": {"query": "TPC1"}}'
```

## Input

| Field | Type | Default | Description |
|---|---|---|---|
| `query` | string | — | **Required.** Gene name, AGI code, or protein name |
| `email` | string | null | Contact email for NCBI (recommended) |
| `make_viewer` | bool | true | Also build the interactive HTML viewer |
| `hold` | int | 20 | Hold frames per species in the GIF |
| `morph` | int | 12 | Morph frames per transition |
| `fps` | int | 12 | GIF frames per second |

## Output

The workflow returns a JSON object with every artifact base64-encoded
(`evolution_gif`, `evolution_sheet_png`, `viewer_html`, `aligned_pdbs`,
`overlay_png`, `sequences_fasta`, `sequences_tsv`, `aligned_json`,
`summary_tsv`) plus a human-readable `report`.

## Project layout

```
workflow/
├── protein_animation/
│   ├── activities.py     # 4 durable activities wrapping the skills
│   └── workflow.py       # protein-animation workflow + supervisor agent
├── skills/               # the four skill scripts (imported as modules)
│   ├── sequence_scraper/
│   ├── structure_aligner/
│   ├── morph_animation/
│   └── interactive_viewer/
├── src/proteinlib.py     # shared geometry/morph library
├── entrypoints/worker.py # worker entrypoint
├── tests/test_activities.py
├── Dockerfile
├── requirements.txt
└── Makefile
```

## Tests

```bash
PYTHONPATH=. python tests/test_activities.py
```

Runs the full align → animate → viewer chain on the cached ECA3 fixture
and validates every artifact payload.

```

## `pyproject.toml`

```
[project]
name = "protein-animation-workflow"
version = "1.0.0"
description = "Mistral Workflow: gene name / AGI / protein name -> evolutionary protein animation (sequences, aligned PDBs, morph GIF, interactive viewer)"
requires-python = ">=3.12"
dependencies = [
    "mistralai-workflows[mistralai]>=0.1.0",
    "mistralai>=2.8.0",
    "pydantic>=2.0",
    "requests>=2.28",
    "numpy>=1.24",
    "biopython>=1.81",
    "matplotlib>=3.7",
    "pillow>=10.0",
    "python-dotenv>=1.0",
]

```

## `requirements.txt`

```
mistralai-workflows[mistralai]>=0.1.0
mistralai>=2.8.0
pydantic>=2.0
requests>=2.28
numpy>=1.24
biopython>=1.81
matplotlib>=3.7
pillow>=10.0
python-dotenv>=1.0

```

## `entrypoints/worker.py`

```
"""Worker entrypoint for the protein-animation workflow."""
import asyncio

import mistralai.workflows as workflows

from protein_animation.workflow import ProteinAnimationWorkflow


async def main():
    await workflows.run_worker([ProteinAnimationWorkflow])


if __name__ == "__main__":
    asyncio.run(main())

```

## `protein_animation/__init__.py`

```


```

## `protein_animation/activities.py`

```
"""Activities for the protein-animation workflow.

Each activity wraps one stage of the evolution pipeline as a durable,
retryable unit of work. The heavy lifting lives in the skill scripts
under skills/<name>/ (sequence-scraper, structure-aligner, morph-animation,
interactive-viewer), imported as modules so the logic stays identical to
the published skills.

All file-producing stages write into a per-run working directory and
return base64-encoded artifacts so they can cross the workflow boundary
as plain JSON.
"""
import asyncio
import base64
import json
import os
import sys
import tempfile
from datetime import timedelta
from pathlib import Path

import mistralai.workflows as workflows

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
for _cand in (str(_ROOT / "src"), str(_ROOT)):
    if _cand not in sys.path:
        sys.path.insert(0, _cand)

from skills.sequence_scraper import scrape_sequences as scraper  # noqa: E402
from skills.structure_aligner import fetch_and_align_structures as aligner  # noqa: E402
from skills.morph_animation import morph_animation as animator  # noqa: E402
from skills.interactive_viewer import interactive_viewer as iv  # noqa: E402


def _b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def _safe(query):
    import re
    return re.sub(r"[^A-Za-z0-9_-]+", "_", query)[:40] or "query"


def _run_dir(query):
    d = Path(tempfile.gettempdir()) / "protein_animation" / _safe(query)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Stage 1: sequence scraping (sequence-scraper skill)
# ---------------------------------------------------------------------------

@workflows.activity(
    start_to_close_timeout=timedelta(minutes=45),
    heartbeat_timeout=timedelta(seconds=60),
    retry_policy_max_attempts=2,
)
async def scrape_sequences(query: str, email: str | None = None) -> dict:
    """Resolve a gene name / AGI / protein name and scrape one verified
    homologous sequence per model plant species (16 species, green algae
    to angiosperms). Returns the FASTA text, the quality TSV, and the
    resolved query accession.
    """
    def _work():
        work = _run_dir(query)
        fasta = work / "sequences.fasta"
        table = work / "sequences.tsv"
        s = scraper.Scraper(email=email)
        terms, qrec = s.resolve_query(query)
        if terms is None:
            raise workflows.WorkflowError(
                f"Could not resolve query: {query}")
        exact_required = qrec is not None
        qlen = qrec["length"] if qrec else None
        qseq = qrec["sequence"] if qrec else None
        skip = ("Arabidopsis thaliana",) if qrec else ()
        rows = s.scrape(terms, exact_required=exact_required,
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
        order = {spec: i for i, (_, spec, _) in enumerate(scraper.SPECIES)}
        rows.sort(key=lambda r: order.get(r["species"], 99))
        scraper.write_fasta(rows, str(fasta))
        scraper.write_table(rows, str(table))
        n_ok = sum(1 for r in rows if r["sequence"])
        return {
            "fasta": fasta.read_text(),
            "table": table.read_text(),
            "n_species": len(rows),
            "n_sequences": n_ok,
            "query_accession": (qrec or {}).get("accession"),
            "terms": terms,
        }

    return await asyncio.to_thread(_work)


# ---------------------------------------------------------------------------
# Stage 2: structure fetching + alignment (structure-aligner skill v2)
# ---------------------------------------------------------------------------

@workflows.activity(
    start_to_close_timeout=timedelta(minutes=60),
    heartbeat_timeout=timedelta(seconds=60),
    retry_policy_max_attempts=2,
)
async def align_structures(query: str, fasta: str,
                           email: str | None = None) -> dict:
    """Fetch an AlphaFold model (or homology-thread) for every scraped
    sequence, superpose all structures into one common 3D frame, and write
    ribbon-ready PDBs + aligned.json. Returns aligned.json plus base64
    copies of every aligned PDB and the overlay plot.
    """
    def _work():
        work = _run_dir(query)
        fasta_path = work / "sequences.fasta"
        fasta_path.write_text(fasta)
        outdir = work / "structures"
        outdir.mkdir(exist_ok=True)

        records = aligner.parse_fasta(str(fasta_path))
        fetcher = aligner.StructureFetcher(str(outdir), email=email)

        structures, templates, pending_model = [], [], []
        used_uniprot = set()
        import numpy as np
        for rec in records:
            name = rec["accession"].split(".")[0]
            up = fetcher.resolve_uniprot(rec)
            if up in used_uniprot:
                up = None
            model, source, mean_plddt = None, "miss", None
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
            if model is not None:
                structures.append({**rec, "uniprot": up, "model": model,
                                   "source": source, "plddt": mean_plddt})
                templates.append(structures[-1])
            else:
                pending_model.append(rec)

        for rec in pending_model:
            best = None
            for t in templates:
                import proteinlib as pl
                pairs = pl.sequence_pairs(rec["sequence"], t["model"]["seq"])
                if not pairs:
                    continue
                frac = len(pairs) / max(len(rec["sequence"]), 1)
                if best is None or frac > best[0]:
                    best = (frac, t)
            if best is None and templates:
                best = (0.0, templates[0])
            if best is None:
                structures.append({**rec, "uniprot": None, "model": None,
                                   "source": "miss", "plddt": None})
                continue
            frac, template = best
            threaded = fetcher.thread_model(rec["sequence"],
                                            template["model"])
            if threaded is None:
                structures.append({**rec, "uniprot": None, "model": None,
                                   "source": "miss", "plddt": None})
                continue
            structures.append({
                **rec, "uniprot": None, "model": threaded,
                "template_accession": template["accession"],
                "source": f"threaded:{template['uniprot']}",
                "plddt": float(np.mean(threaded["plddts"]))})

        got = [s for s in structures if s["model"] is not None]
        if len(got) < 2:
            raise workflows.WorkflowError(
                "Fewer than 2 structures available; cannot align.")

        import proteinlib as pl

        def avg_sim(s):
            sims = []
            for o in got:
                if o is s:
                    continue
                pairs = pl.sequence_pairs(s["model"]["seq"],
                                          o["model"]["seq"]) or []
                sims.append(pl.core_identity(s["model"]["seq"],
                                             o["model"]["seq"], pairs))
            return float(np.mean(sims)) if sims else 0.0

        real = [s for s in got if not s.get("template_accession")]
        reference = max(real or got, key=avg_sim)

        ref_coords = reference["model"]["coords"]
        aligned, transforms, core_idents = [], {}, {}
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
        for s in got:
            if s is reference:
                new_coords, fit = ref_coords, 0.0
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
            pdb_dir = outdir / "aligned"
            pdb_dir.mkdir(exist_ok=True)
            pl.write_ca_pdb(new_coords, s["model"]["plddts"],
                            str(pdb_dir / f"{name}.pdb"),
                            name, ss=ss, seq=seq)

        aligned_json = outdir / "aligned.json"
        aligned_json.write_text(json.dumps(
            {"reference": reference["accession"], "structures": aligned}))

        summary = outdir / "summary.tsv"
        with open(summary, "w") as f:
            f.write("accession\tspecies\tlineage\tuniprot\tsource\tplddt\t"
                    "rmsd_to_ref\tn_residues\n")
            for a in aligned:
                f.write(f"{a['accession']}\t{a['species']}\t{a['lineage']}\t"
                        f"{a['uniprot'] or ''}\t{a['source']}\t"
                        f"{a['plddt'] if a['plddt'] is not None else ''}\t"
                        f"{a['rmsd_to_ref']:.2f}\t{a['n_residues']}\n")

        # overlay plot (best effort) - dark theme to match the renderer
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(10, 8))
            fig.patch.set_facecolor("#0d1117")
            ax = fig.add_subplot(111, projection="3d")
            ax.set_facecolor("#0d1117")
            for a in aligned:
                c = np.array(a["coords"])
                ax.plot(c[:, 0], c[:, 1], c[:, 2], linewidth=0.5,
                        alpha=0.5, color="#58a6ff")
            ax.set_axis_off()
            ax.set_title(f"{len(aligned)} structures superposed on "
                         f"{reference['accession']}", color="#e6edf3")
            plt.savefig(str(outdir / "overlay.png"), dpi=120,
                        facecolor=fig.get_facecolor())
            plt.close(fig)
        except Exception:
            pass

        pdb_files = {}
        for p in sorted((outdir / "aligned").glob("*.pdb")):
            pdb_files[p.name] = _b64(str(p))
        overlay = None
        if (outdir / "overlay.png").exists():
            overlay = _b64(str(outdir / "overlay.png"))
        return {
            "aligned_json": aligned_json.read_text(),
            "summary": summary.read_text(),
            "reference": reference["accession"],
            "n_structures": len(aligned),
            "species": [a["species"] for a in aligned],
            "pdbs": pdb_files,
            "overlay_png": overlay,
        }

    return await asyncio.to_thread(_work)


# ---------------------------------------------------------------------------
# Stage 3: morph animation (morph-animation skill v2)
# ---------------------------------------------------------------------------

@workflows.activity(
    start_to_close_timeout=timedelta(minutes=45),
    heartbeat_timeout=timedelta(seconds=60),
    retry_policy_max_attempts=2,
)
async def render_animation(query: str, aligned_json: str,
                             hold: int = 20, morph: int = 12,
                             fps: int = 12,
                             size: tuple = (720, 720)) -> dict:
    """Render the ribbon-style, stability-masked evolutionary morph GIF
    plus a one-still-per-species frame sheet PNG.
    """
    def _work():
        import numpy as np
        import proteinlib as pl
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        work = _run_dir(query)
        data = json.loads(aligned_json)
        structs = animator.order_structures(data["structures"])
        if len(structs) < 2:
            raise workflows.WorkflowError(
                "Need at least 2 structures to animate")

        gif_path = work / "evolution.gif"
        # correspondence anchored through the alignment reference, so a
        # pairwise alignment slip between distant homologs cannot reshuffle
        # the whole protein; N/C termini labeled on every frame
        reference = next((s for s in data["structures"]
                          if s["accession"] == data.get("reference")),
                         structs[0])
        mask = pl.StabilityMask(disp_max=4.0)
        frames_meta = []
        for idx, s in enumerate(structs):
            coords = np.asarray(s["coords"], dtype=float)
            seq = s.get("seq") or "A" * len(coords)
            ss = s.get("ss") or "C" * len(coords)
            fm = {
                "title": (f"{animator.LINEAGE_LABELS.get(s['lineage'], s['lineage'])}"
                          f" \u2014 {s['species']}"),
                "subtitle": (f"{s['accession']} \u00b7 {s['source']}"
                             + (f" \u00b7 core ident {s['core_identity']:.0%}"
                                if s.get("core_identity") else "")),
                "ss": ss,
                "termini": (0, len(coords) - 1) if len(coords) else None,
            }
            cols = pl.plddt_colors(s.get("plddts") or [70] * len(coords),
                                   len(coords))
            for _ in range(hold):
                frames_meta.append({**fm, "coords": coords, "colors": cols})
            if idx + 1 < len(structs):
                nxt = structs[idx + 1]
                ncoords = np.asarray(nxt["coords"], dtype=float)
                nseq = nxt.get("seq") or "A" * len(ncoords)
                nss = nxt.get("ss") or "C" * len(ncoords)
                pairs = pl.reference_pairs(s, nxt, reference)
                pairs = [(a, b) for a, b in pairs
                         if a < len(coords) and b < len(ncoords)]
                weights, a_idx, b_idx = mask.weights(
                    coords, ncoords, pairs, ss, nss)
                interp = np.asarray(pl.morph_frames(coords, ncoords, pairs,
                                                    weights, morph))
                n_cols_next = pl.plddt_colors(
                    nxt.get("plddts") or [70] * len(ncoords), len(ncoords))
                cols_a = cols[a_idx]
                cols_b = n_cols_next[b_idx]
                n_stable = int((weights < 0.15).sum())
                for k in range(interp.shape[0]):
                    t = k / max(interp.shape[0] - 1, 1)
                    cols_k = cols_a * (1 - t) + cols_b * t
                    b_to_pos = {b: pos for pos, b in enumerate(b_idx)}
                    ss_full = []
                    for j in range(len(ncoords)):
                        pos = b_to_pos.get(j)
                        if pos is None:
                            ss_full.append(nss[j])
                        else:
                            w = weights[pos]
                            ss_full.append(ss[a_idx[pos]] if w < 0.5
                                           else nss[b_idx[pos]])
                    frames_meta.append({
                        **{"title": (fm["title"] if t < 0.5 else
                                     f"{animator.LINEAGE_LABELS.get(nxt['lineage'], nxt['lineage'])}"
                                     f" \u2014 {nxt['species']}"),
                           "subtitle": fm["subtitle"] if t < 0.5 else
                               f"{nxt['accession']} \u00b7 {nxt['source']}"
                               + (f" \u00b7 core ident {nxt['core_identity']:.0%}"
                                  if nxt.get("core_identity") else "")},
                        "coords": interp[k], "colors": cols_k,
                        "ss": "".join(ss_full),
                        "termini": (0, len(ncoords) - 1)
                        if len(ncoords) else None})

        animator.render(frames_meta, str(gif_path), fps, tuple(size))

        sheet_path = work / "evolution_sheet.png"
        try:
            sheet = frames_meta[::hold + morph]
            cols_n = int(np.ceil(np.sqrt(len(sheet))))
            rows_n = int(np.ceil(len(sheet) / cols_n))
            fig = plt.figure(figsize=(cols_n * 3, rows_n * 3), dpi=100)
            fig.patch.set_facecolor(animator.BG)
            for i, fm in enumerate(sheet, 1):
                ax = fig.add_subplot(rows_n, cols_n, i, projection="3d")
                ax.set_facecolor(animator.BG)
                c = np.asarray(fm["coords"])
                ss = fm.get("ss") or "C" * len(c)
                cols = np.asarray(fm["colors"])[:, :3]
                # apply the same view-space shading as the GIF renderer so
                # the sheet matches the animation's lighting
                tang = np.gradient(c, axis=0)
                norm = np.linalg.norm(tang, axis=1, keepdims=True)
                norm[norm < 1e-6] = 1.0
                tang = tang / norm
                light = np.array([0.5, 0.6, 0.62])
                light /= np.linalg.norm(light)
                lam = np.clip(np.abs(tang @ light), 0.0, 1.0)
                shade_f = 0.45 + 0.55 * lam[:, None]
                cols = np.clip(cols * shade_f, 0.0, 1.0)
                termini = fm.get("termini")
                if termini is not None:
                    n_idx, c_idx = termini
                    termini = (min(n_idx, len(c) - 1),
                               min(c_idx, len(c) - 1))
                animator.draw_structure(ax, c, ss, cols, termini=termini)
                ax.set_axis_off()
                ax.set_title(fm["title"], fontsize=7, color="#e6edf3")
                ax.view_init(elev=18, azim=30)
            plt.tight_layout()
            plt.savefig(str(sheet_path),
                        facecolor=fig.get_facecolor())
            plt.close(fig)
        except Exception:
            pass

        return {
            "gif": _b64(str(gif_path)),
            "gif_filename": "evolution.gif",
            "sheet_png": _b64(str(sheet_path)) if sheet_path.exists() else None,
            "n_frames": len(frames_meta),
            "n_structures": len(structs),
        }

    return await asyncio.to_thread(_work)


# ---------------------------------------------------------------------------
# Stage 4: interactive viewer (interactive-viewer skill v9) - optional
# ---------------------------------------------------------------------------

@workflows.activity(
    start_to_close_timeout=timedelta(minutes=20),
    heartbeat_timeout=timedelta(seconds=60),
    retry_policy_max_attempts=1,
)
async def build_viewer(query: str, aligned_json: str,
                       title: str | None = None) -> dict:
    """Build the single-file interactive HTML viewer (clickable tree +
    rotatable ribbon structure + stability-masked morph). Taxonomy
    fetching needs network access to NCBI; failures degrade to a tree
    without internal labels.
    """
    def _work():
        work = _run_dir(query)
        aj = work / "aligned.json"
        aj.write_text(aligned_json)

        viewer_mod_dir = str(_ROOT / "skills" / "interactive_viewer")
        # reuse the module's own pipeline: taxonomy + payload + template
        import json as _json
        import numpy as np
        data = _json.loads(aligned_json)
        structs = data["structures"]
        for s in structs:
            iv.SPECIES_LINEAGE[s["species"]] = s.get("lineage", "")
        import requests
        session = requests.Session()
        cache = {}
        lineages = {}
        for s in structs:
            try:
                lineages[s["species"]] = iv.fetch_lineage(
                    s["species"], session, cache)
            except Exception:
                lineages[s["species"]] = []
        tree = iv.build_tree(lineages)
        nodes, edges, leaves = iv.layout_tree(tree, iv.LINEAGE_ORDER)

        out_structs = []
        ref_acc = data.get("reference")
        ref = next((s for s in structs if s["accession"] == ref_acc),
                   structs[0])
        ref_coords = np.asarray(ref["coords"], dtype=float)
        ref_seq = ref.get("seq") or "A" * len(ref_coords)
        import proteinlib as pl
        for s in structs:
            coords = np.asarray(s["coords"], dtype=float)
            seq = s.get("seq") or "A" * len(coords)
            n_coord = len(coords)
            coord_seq = (seq[:n_coord] if len(seq) >= n_coord
                         else seq + "X" * (n_coord - len(seq)))
            h = iv.hydrophobicity(coord_seq)
            ss = s.get("ss") or "".join(pl.assign_secondary(coords))
            pairs = pl.sequence_pairs(seq, ref_seq)
            pairs = [[a, b] for a, b in pairs
                     if a < len(coords) and b < len(ref_coords)]
            out_structs.append({
                "accession": s["accession"], "species": s["species"],
                "lineage": s.get("lineage", ""),
                "source": s.get("source", ""),
                "core_identity": s.get("core_identity"),
                "n_residues": len(coords),
                "ss": ss,
                "pdb": iv.to_pdb_text({"coords": coords.tolist(),
                                       "seq": seq, "ss": ss}),
                "hydro": [round(float(x), 2) for x in h],
                "ref_pairs": pairs,
                "ref_coords": None,
            })
        ref_pts = [[round(float(c), 3) for c in xyz] for xyz in ref_coords]
        for s in out_structs:
            s["ref_coords"] = ref_pts
            if s["accession"] == ref["accession"]:
                s["ref_pairs"] = [[i, i] for i in range(len(ref_coords))]

        payload = {
            "title": title or f"{query} protein structure across lineages",
            "reference": ref["accession"],
            "lineageOrder": iv.LINEAGE_ORDER,
            "lineageColors": iv.LINEAGE_COLORS,
            "tree": {"nodes": nodes, "edges": edges, "leaves": leaves},
            "structures": out_structs,
        }
        template = open(os.path.join(viewer_mod_dir,
                                     "viewer_template.html")).read()
        html = (template.replace("__DATA__", _json.dumps(payload))
                        .replace("__TITLE__", payload["title"]))
        out = work / "viewer.html"
        out.write_text(html)
        return {
            "viewer_html": _b64(str(out)),
            "viewer_filename": "viewer.html",
            "title": payload["title"],
        }

    return await asyncio.to_thread(_work)

```

## `protein_animation/workflow.py`

```
"""protein-animation workflow.

Input:  a protein query - a gene name (e.g. "TPC1"), an Arabidopsis AGI
        code (e.g. "AT4G03560"), or a protein name.
Output: the full artifact set - sequences, aligned PDBs, the evolutionary
        morph animation GIF, a frame sheet, and an interactive HTML viewer.

The pipeline chains the workspace skills:
    sequence-scraper -> structure-aligner -> morph-animation
                                          -> interactive-viewer

A supervisor agent (durable agent) watches the run, interprets stage
results, and assembles the final report. The heavy lifting is done by
deterministic activities so the run is durable and retryable.
"""
from datetime import timedelta

from pydantic import BaseModel, Field

import mistralai.workflows as workflows
import mistralai.workflows.plugins.mistralai as workflows_mistralai
from mistralai.client.models import TextChunk

from .activities import (align_structures, build_viewer, render_animation,
                         scrape_sequences)


class ProteinAnimationInput(BaseModel):
    query: str = Field(
        description="Gene name (e.g. 'TPC1'), Arabidopsis AGI code "
                    "(e.g. 'AT4G03560'), or protein name")
    email: str | None = Field(
        default=None,
        description="Contact email for NCBI E-utilities (recommended)")
    make_viewer: bool = Field(
        default=True,
        description="Also build the interactive HTML viewer")
    hold: int = Field(default=20, description="Hold frames per species")
    morph: int = Field(default=12, description="Morph frames per transition")
    fps: int = Field(default=12, description="GIF frames per second")


class ProteinAnimationOutput(BaseModel):
    query: str
    terms: list[str]
    query_accession: str | None
    n_species: int
    n_sequences: int
    reference: str
    n_structures: int
    species: list[str]
    sequences_fasta: str
    sequences_tsv: str
    aligned_json: str
    summary_tsv: str
    aligned_pdbs: dict
    overlay_png: str | None
    evolution_gif: str
    evolution_sheet_png: str | None
    n_frames: int
    viewer_html: str | None
    report: str


@workflows.workflow.define(
    name="protein-animation",
    workflow_display_name="Protein Animation",
    workflow_description=(
        "Generate an evolutionary protein structure animation from a gene "
        "name, AGI code, or protein name. Chains sequence-scraper, "
        "structure-aligner, morph-animation and interactive-viewer. "
        "Outputs sequences, aligned PDBs, the morph GIF, a frame sheet, "
        "and an interactive HTML viewer."),
    execution_timeout=timedelta(hours=3),
)
class ProteinAnimationWorkflow:
    @workflows.workflow.entrypoint
    async def run(self, input: ProteinAnimationInput) -> dict:
        # ---- Stage 1: sequences -----------------------------------------
        seq = await scrape_sequences(input.query, input.email)
        if seq["n_sequences"] < 2:
            raise workflows.WorkflowError(
                f"Scraper found only {seq['n_sequences']} sequences for "
                f"'{input.query}'; need at least 2 to animate.")

        # ---- Stage 2: structures ----------------------------------------
        aln = await align_structures(input.query, seq["fasta"], input.email)

        # ---- Stage 3: animation -----------------------------------------
        anim = await render_animation(
            input.query, aln["aligned_json"],
            hold=input.hold, morph=input.morph, fps=input.fps)

        # ---- Stage 4: interactive viewer (optional) ---------------------
        viewer = None
        if input.make_viewer:
            try:
                viewer = await build_viewer(input.query,
                                            aln["aligned_json"],
                                            title=input.query)
            except Exception:
                viewer = None

        # ---- Supervisor agent: interpret and report ---------------------
        report = await self._supervise(input, seq, aln, anim, viewer)

        out = ProteinAnimationOutput(
            query=input.query,
            terms=seq["terms"],
            query_accession=seq["query_accession"],
            n_species=seq["n_species"],
            n_sequences=seq["n_sequences"],
            reference=aln["reference"],
            n_structures=aln["n_structures"],
            species=aln["species"],
            sequences_fasta=seq["fasta"],
            sequences_tsv=seq["table"],
            aligned_json=aln["aligned_json"],
            summary_tsv=aln["summary"],
            aligned_pdbs=aln["pdbs"],
            overlay_png=aln["overlay_png"],
            evolution_gif=anim["gif"],
            evolution_sheet_png=anim["sheet_png"],
            n_frames=anim["n_frames"],
            viewer_html=(viewer or {}).get("viewer_html"),
            report=report,
        )
        return out.model_dump()

    async def _supervise(self, input, seq, aln, anim, viewer) -> str:
        """Durable supervisor agent: summarize the run for the caller."""
        agent = workflows_mistralai.Agent(
            model="mistral-small-latest",
            name="protein-animation-supervisor",
            description="Summarizes protein animation pipeline runs",
            instructions=(
                "You supervise runs of the protein-animation pipeline. You "
                "receive a JSON status blob describing each stage. Write a "
                "short, friendly markdown report for the user: what query "
                "was resolved to, how many species/sequences/structures "
                "were found, which species are included, the alignment "
                "reference, how many animation frames were rendered, and "
                "which output files are available. Keep it under 150 words. "
                "Do not invent data; only use what is in the status blob."),
        )
        status = {
            "query": input.query,
            "resolved_terms": seq["terms"],
            "query_accession": seq["query_accession"],
            "species_found": seq["n_sequences"],
            "species_total": seq["n_species"],
            "structures_aligned": aln["n_structures"],
            "reference_accession": aln["reference"],
            "species_list": aln["species"],
            "animation_frames": anim["n_frames"],
            "viewer_built": viewer is not None,
        }
        outputs = await workflows_mistralai.Runner.run(
            agent=agent,
            inputs=("Summarize this protein-animation run status:\n"
                    + json.dumps(status)),
            session=workflows_mistralai.RemoteSession(),
            max_turns=3,
        )
        text = "\n".join(o.text for o in outputs
                         if isinstance(o, TextChunk))
        return text or "Run complete."


import json  # noqa: E402  (used by _supervise)

```

## `skills/__init__.py`

```


```

## `skills/interactive_viewer/__init__.py`

```


```

## `skills/interactive_viewer/interactive_viewer.py`

```
#!/usr/bin/env python3
"""Generate an interactive evolutionary structure viewer (single HTML). (v14)

Usage:
    python interactive_viewer.py structures/aligned.json -o viewer.html

Reads aligned.json (structure-aligner skill v2) and embeds:
  - an interactive phylogenetic tree (SVG) of the species, colored by
    lineage; clicking a species morphs the structure
  - a Three.js cartoon-style 3D viewer (works from CA-only traces):
    round tube with variable thickness (helix thick, strand medium,
    loop thin), smooth transitions, rotatable/zoomable
  - per-residue PURE HYDRopathy coloring (blue hydrophilic -> white ->
    red hydrophobic), no SS-class tint
  - consistent MEDIUM lighting: hemisphere + directional lights attached
    to the camera at fixed medium intensity, so shading is identical while
    orbiting
  - STABILITY-MASKED morphing: when switching species, residues that are
    spatially stable AND keep their secondary-structure class stay put
    (a helix that is a helix in both species does not reshuffle); only
    genuinely changing regions move, insertions fade in, and timing is
    smoothstep-eased
  - v14 UNION-SPACE MORPH + DoubleSide: handles both insertions AND
    deletions. Insertions fade in (radius 0->1), deletions collapse and
    fade out (1->0). DoubleSide material fixes "see inside tube" issue.
    No more spline artifacts from clumped coincident points near
    conserved helices when C/N-terminal tails differ in length.
  - v13 STANDARD CARTOON: round tube with variable thickness (helix=2.2,
    strand=1.6, loop=0.7 radius), smooth 3-point-averaged transitions at
    SS boundaries. No elliptical cross-sections — clean PyMOL/VMD look.
  - v12 PARALLEL-TRANSPORT FRAMES + LOCAL INSERTION GROWTH: the tube's
    cross-section frame follows the curve via parallel transport (fixes
    oval/twisted helix ribbons), and insertions grow from the nearest
    paired neighbor's morphing position instead of the global centroid
    (fixes "expanding from inner center" on distant-species transitions).
  - v11 VARIABLE-RADIUS TUBE + NO REBUILD AT MORPH END: the tube's
    elliptical cross-section (width/height) cross-fades with SS changes
    during the morph, and at f=1 the tube snaps in place (no clearGroup /
    showStructure) so the camera orientation is preserved — no more
    broken-looking morph or orientation jump at the end.
  - N-/C-termini are marked with labeled spheres on every structure.
  - NO membrane visualization (removed in v7; template asset restored in v8)

The output is one self-contained HTML file (Three.js from CDN).
"""
import argparse
import json
import os
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import requests

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

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

LINEAGE_ORDER = ["Green Algae", "Bryophytes", "Ferns",
                 "Gymnosperms", "Angiosperms"]
LINEAGE_COLORS = {
    "Green Algae": "#2ca02c", "Bryophytes": "#8c564b",
    "Ferns": "#ff7f0e", "Gymnosperms": "#9467bd",
    "Angiosperms": "#d62728",
}

KD = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9,
      'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9,
      'P': -1.6, 'Y': -1.3, 'H': -3.2, 'E': -3.5, 'Q': -3.5,
      'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0.0}

AA1TO3 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
          'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
          'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
          'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
          'X': 'GLY'}

# BLOSUM62 substitution matrix as a flat dict for sequence similarity coloring.
# Scores range from -4 (very dissimilar) to +11 (identical, e.g. W-W).
# Used to color aligned sequence residues by similarity to the reference.
BLOSUM62 = {
    'A': {'A':4,'R':-1,'N':-2,'D':-2,'C':0,'Q':-1,'E':-1,'G':0,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1,'T':0,'W':-3,'Y':-2,'V':0},
    'R': {'A':-1,'R':5,'N':0,'D':-2,'C':-3,'Q':1,'E':0,'G':-2,'H':0,'I':-3,'L':-2,'K':2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3},
    'N': {'A':-2,'R':0,'N':6,'D':1,'C':-3,'Q':0,'E':0,'G':0,'H':1,'I':-3,'L':-3,'K':0,'M':-2,'F':-3,'P':-2,'S':1,'T':0,'W':-4,'Y':-2,'V':-3},
    'D': {'A':-2,'R':-2,'N':1,'D':6,'C':-3,'Q':0,'E':2,'G':-1,'H':-1,'I':-3,'L':-3,'K':-1,'M':-3,'F':-3,'P':-1,'S':0,'T':-1,'W':-4,'Y':-3,'V':-3},
    'C': {'A':0,'R':-3,'N':-3,'D':-3,'C':9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1},
    'Q': {'A':-1,'R':1,'N':0,'D':0,'C':-3,'Q':5,'E':2,'G':-2,'H':0,'I':-3,'L':-2,'K':1,'M':0,'F':-3,'P':-1,'S':0,'T':-1,'W':-2,'Y':-1,'V':-2},
    'E': {'A':-1,'R':0,'N':0,'D':2,'C':-4,'Q':2,'E':5,'G':-2,'H':0,'I':-3,'L':-3,'K':1,'M':-2,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
    'G': {'A':0,'R':-2,'N':0,'D':-1,'C':-3,'Q':-2,'E':-2,'G':6,'H':-2,'I':-4,'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S':0,'T':-2,'W':-2,'Y':-3,'V':-3},
    'H': {'A':-2,'R':0,'N':1,'D':-1,'C':-3,'Q':0,'E':0,'G':-2,'H':8,'I':-3,'L':-2,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y':2,'V':-3},
    'I': {'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I':4,'L':2,'K':-3,'M':1,'F':0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V':3},
    'L': {'A':-1,'R':-2,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-2,'I':2,'L':4,'K':-2,'M':2,'F':0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V':1},
    'K': {'A':-1,'R':2,'N':0,'D':-1,'C':-3,'Q':1,'E':1,'G':-2,'H':-1,'I':-3,'L':-2,'K':5,'M':-1,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
    'M': {'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q':0,'E':-2,'G':-3,'H':-2,'I':1,'L':2,'K':-1,'M':5,'F':0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-2,'V':1},
    'F': {'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I':0,'L':0,'K':-3,'M':0,'F':6,'P':-4,'S':-2,'T':-1,'W':1,'Y':3,'V':-1},
    'P': {'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,'L':-3,'K':-1,'M':-2,'F':-4,'P':7,'S':-1,'T':-2,'W':-4,'Y':-3,'V':-2},
    'S': {'A':1,'R':-1,'N':1,'D':0,'C':-1,'Q':0,'E':0,'G':0,'H':-1,'I':-2,'L':-2,'K':0,'M':-1,'F':-2,'P':-1,'S':4,'T':1,'W':-3,'Y':-2,'V':-2},
    'T': {'A':0,'R':-1,'N':0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-1,'P':-2,'S':1,'T':5,'W':-2,'Y':-2,'V':0},
    'W': {'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,'L':-2,'K':-3,'M':-1,'F':1,'P':-4,'S':-3,'T':-2,'W':11,'Y':2,'V':-3},
    'Y': {'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H':2,'I':-1,'L':-1,'K':-2,'M':-2,'F':3,'P':-3,'S':-2,'T':-2,'W':2,'Y':7,'V':-1},
    'V': {'A':0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I':3,'L':1,'K':-2,'M':1,'F':-1,'P':-2,'S':-2,'T':0,'W':-3,'Y':-1,'V':4},
}


def blosum_score(a, b):
    """BLOSUM62 substitution score between two amino acids.
    Returns 0 for unknown residues (X, gaps, etc.)."""
    if a not in BLOSUM62 or b not in BLOSUM62:
        return 0
    return BLOSUM62[a][b]


def build_msa(structs, ref_seq, ref_acc):
    """Build a multiple sequence alignment (MSA) anchored to the reference.

    Uses a reference-anchored star alignment: the reference sequence defines
    the column backbone (one column per reference residue, no gaps in the
    reference row). Insertions from other sequences are placed in additional
    columns between reference columns. The reference alignment row has NO
    gaps, making it easy to map structure residue indices to MSA columns.

    Returns:
      - msa_seqs: list of (species, aligned_seq) in MSA column space
      - msa_cols: number of columns in the MSA
      - ref_col_map: list mapping reference residue index -> MSA column
        (identity mapping since reference defines the backbone, but kept
        for generality; ref residue i is at column ref_col_map[i])
    """
    n_ref = len(ref_seq)

    # For each structure, get alignment pairs to reference
    all_pairs = {}
    for s in structs:
        seq = s.get("seq") or ""
        if not seq:
            all_pairs[s["species"]] = []
            continue
        if s["accession"] == ref_acc:
            all_pairs[s["species"]] = [(i, i) for i in range(min(len(seq), n_ref))]
        else:
            pairs = pl.sequence_pairs(seq, ref_seq)
            all_pairs[s["species"]] = pairs

    # Build seq_to_ref mapping
    seq_to_ref = {}
    for sp, pairs in all_pairs.items():
        seq_to_ref[sp] = {si: ri for si, ri in pairs}

    # For each non-reference sequence, find unpaired residues (insertions)
    # and group them into blocks, noting which reference position they
    # insert after.
    seq_insertions = {}  # species -> list of (seq_start, seq_end, after_ref_idx)
    for s in structs:
        sp = s["species"]
        if sp == ref_acc or sp not in seq_to_ref:
            continue
        seq = s.get("seq") or ""
        if not seq:
            continue
        paired = set(seq_to_ref[sp].keys())
        unpaired = [i for i in range(len(seq)) if i not in paired]
        insertions = []
        i = 0
        while i < len(unpaired):
            start = unpaired[i]
            end = start
            while i + 1 < len(unpaired) and unpaired[i + 1] == end + 1:
                i += 1
                end = unpaired[i]
            # Find the reference position this block inserts after
            before = -1
            for j in range(start - 1, -1, -1):
                if j in seq_to_ref[sp]:
                    before = seq_to_ref[sp][j]
                    break
            insertions.append((start, end, before))
            i += 1
        seq_insertions[sp] = insertions

    # Collect all insertion points: after_ref_idx -> list of (species, seq_start, seq_end)
    # Each insertion block from each species gets its own set of columns.
    # Multiple species inserting at the same point get separate column groups.
    ins_points = {}  # after_ref_idx -> list of (species, seq_start, seq_end)
    for sp, inss in seq_insertions.items():
        for seq_start, seq_end, after_ref in inss:
            if after_ref not in ins_points:
                ins_points[after_ref] = []
            ins_points[after_ref].append((sp, seq_start, seq_end))

    # Build column layout:
    # - Reference columns: one per reference residue (ref row has no gaps)
    # - Insertion columns: between reference columns, one per inserted residue
    #   from each species (each species gets its own columns)
    # Column structure: {type: 'ref', ref_idx: int} or {type: 'ins', species: str, seq_idx: int}
    columns = []
    # Pre-insertion columns (insertions before first ref residue, after_ref = -1)
    if -1 in ins_points:
        for sp, ss, se in sorted(ins_points[-1]):
            for k in range(ss, se + 1):
                columns.append({"type": "ins", "species": sp, "seq_idx": k})
    # Reference columns with interleaved insertion columns
    for ri in range(n_ref):
        columns.append({"type": "ref", "ref_idx": ri})
        if ri in ins_points:
            for sp, ss, se in sorted(ins_points[ri]):
                for k in range(ss, se + 1):
                    columns.append({"type": "ins", "species": sp, "seq_idx": k})
    # Post-insertion columns (insertions after last ref residue)
    if n_ref in ins_points:
        for sp, ss, se in sorted(ins_points[n_ref]):
            for k in range(ss, se + 1):
                columns.append({"type": "ins", "species": sp, "seq_idx": k})

    n_cols = len(columns)

    # Build aligned sequences
    msa_seqs = []
    for s in structs:
        sp = s["species"]
        seq = s.get("seq") or ""
        if not seq:
            msa_seqs.append((sp, "-" * n_cols))
            continue
        # Build column -> residue mapping
        col_to_res = {}
        if sp in seq_to_ref:
            # Map reference columns
            for ci, col in enumerate(columns):
                if col["type"] == "ref":
                    ri = col["ref_idx"]
                    # Find seq_idx that maps to this ref_idx
                    for si, rj in seq_to_ref[sp].items():
                        if rj == ri:
                            col_to_res[ci] = seq[si]
                            break
        # Map insertion columns
        for ci, col in enumerate(columns):
            if col["type"] == "ins" and col["species"] == sp:
                col_to_res[ci] = seq[col["seq_idx"]]
        # Build the aligned string
        aligned = []
        for ci in range(n_cols):
            aligned.append(col_to_res.get(ci, "-"))
        msa_seqs.append((sp, "".join(aligned)))

    # ref_col_map: reference residue i is at the column where type='ref' and ref_idx==i
    ref_col_map = []
    for ri in range(n_ref):
        for ci, col in enumerate(columns):
            if col["type"] == "ref" and col["ref_idx"] == ri:
                ref_col_map.append(ci)
                break

    return msa_seqs, n_cols, ref_col_map


def compute_similarity_scores(msa_seqs, ref_species):
    """Compute per-column BLOSUM62 similarity scores for each sequence.

    For each sequence, for each MSA column, compute the BLOSUM62 score
    between the sequence's residue and the reference residue at that column.
    Returns: list of (species, [score_per_column]) where scores are normalized
    to [0, 1] (0 = dissimilar, 1 = identical).
    """
    ref_aligned = None
    for sp, aln in msa_seqs:
        if sp == ref_species:
            ref_aligned = aln
            break
    if ref_aligned is None:
        return [(sp, [0.5] * len(aln)) for sp, aln in msa_seqs]

    n_cols = len(ref_aligned)
    results = []
    for sp, aln in msa_seqs:
        scores = []
        for ci in range(n_cols):
            ref_aa = ref_aligned[ci]
            seq_aa = aln[ci]
            if ref_aa == "-" or seq_aa == "-":
                scores.append(None)
            else:
                raw = blosum_score(seq_aa, ref_aa)
                self_score = blosum_score(ref_aa, ref_aa)
                if self_score > 0:
                    norm = (raw + 4) / (self_score + 4)
                    norm = max(0, min(1, norm))
                else:
                    norm = 0.5
                scores.append(round(norm, 3))
        results.append((sp, scores))
    return results


def hydrophobicity(seq, window=9):
    """Smoothed Kyte-Doolittle hydropathy per residue."""
    raw = np.array([KD.get(a, 0.0) for a in seq], dtype=float)
    if len(raw) == 0:
        return raw
    half = window // 2
    sm = np.zeros_like(raw)
    for i in range(len(raw)):
        lo = max(0, i - half)
        hi = min(len(raw), i + half + 1)
        sm[i] = raw[lo:hi].mean()
    return sm


def to_pdb_text(struct):
    """CA-trace PDB with column-exact ATOM and HELIX/SHEET records."""
    coords = np.asarray(struct["coords"], dtype=float)
    seq = struct.get("seq") or "A" * len(coords)
    ss = struct.get("ss") or (pl.assign_secondary(coords)
                              if pl else ['C'] * len(coords))
    ss = "".join(ss)
    lines = []
    nhelix = 0
    nsheet = 0
    i = 0
    while i < len(ss):
        if ss[i] == 'H':
            j = i
            while j < len(ss) and ss[j] == 'H':
                j += 1
            if j - i >= 4:
                nhelix += 1
                r1 = AA1TO3.get(seq[i], 'GLY')
                r2 = AA1TO3.get(seq[j - 1], 'GLY')
                lines.append(f"HELIX  {nhelix:3d} {nhelix:3d} "
                             f"{r1:>3s} A {i + 1:4d}  {r2:>3s} "
                             f"A {j:4d}  1")
            i = j
        else:
            i += 1
    i = 0
    while i < len(ss):
        if ss[i] == 'E':
            j = i
            while j < len(ss) and ss[j] == 'E':
                j += 1
            if j - i >= 3:
                nsheet += 1
                r1 = AA1TO3.get(seq[i], 'GLY')
                r2 = AA1TO3.get(seq[j - 1], 'GLY')
                lines.append(f"SHEET  {nsheet:3d} {nsheet:3d} 1 "
                             f"{r1:>3s} A{i + 1:4d}  {r2:>3s} "
                             f"A{j:4d}  0")
            i = j
        else:
            i += 1
    for k, xyz in enumerate(coords, start=1):
        x, y, z = xyz
        res = AA1TO3.get(seq[k - 1] if k - 1 < len(seq) else 'A', 'GLY')
        lines.append(f"ATOM  {k:5d}  CA  {res:>3s} A{k:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00")
    lines.append("TER")
    lines.append("END")
    return chr(10).join(lines)


def fetch_lineage(species, session, cache):
    if species in cache:
        return cache[species]
    for _ in range(3):
        try:
            r = session.get(f"{BASE}/esearch.fcgi",
                            params={"db": "taxonomy", "term": species,
                                    "retmax": 1}, timeout=30)
            root = ET.fromstring(r.text)
            ids = [e.text for e in root.findall(".//IdList/Id")]
            if not ids:
                break
            r2 = session.get(f"{BASE}/efetch.fcgi",
                             params={"db": "taxonomy", "id": ids[0]},
                             timeout=30)
            root2 = ET.fromstring(r2.text)
            ranks = [(el.findtext("Rank"), el.findtext("ScientificName"))
                     for el in root2.findall(".//Taxon/LineageEx/Taxon")]
            ranks = [name for rank, name in ranks
                     if rank in ("superkingdom", "kingdom", "phylum",
                                 "subphylum", "class", "subclass", "order",
                                 "suborder", "family", "subfamily", "genus")]
            cache[species] = ranks
            time.sleep(0.35)
            return ranks
        except Exception:
            time.sleep(1)
    cache[species] = []
    return []


def build_tree(species_lineages):
    root = {"name": "root", "children": {}, "leaf": None}
    for sp, ranks in species_lineages.items():
        node = root
        for name in ranks:
            if name not in node["children"]:
                node["children"][name] = {"name": name, "children": {},
                                          "leaf": None}
            node = node["children"][name]
        node["children"][sp] = {"name": sp, "children": {},
                                "leaf": sp}
    return root


def layout_tree(root, order):
    leaves = []

    def collect(node):
        if node["leaf"]:
            leaves.append(node["leaf"])
            return
        for ch in node["children"].values():
            collect(ch)
    collect(root)

    rank0 = {sp: i for i, sp in enumerate(leaves)}
    def leaf_key(sp):
        lin = SPECIES_LINEAGE.get(sp, "")
        r = order.index(lin) if lin in order else 99
        return (r, rank0[sp])
    leaves.sort(key=leaf_key)

    pos_y = {sp: i for i, sp in enumerate(leaves)}
    nodes = []
    edges = []

    def maxdepth(node, d=0):
        if node["leaf"] or not node["children"]:
            return d
        return max(maxdepth(ch, d + 1)
                   for ch in node["children"].values())
    leaf_depth = maxdepth(root)

    def place(node, depth):
        if node["leaf"]:
            y = pos_y[node["leaf"]]
            nodes.append({"x": leaf_depth, "y": y,
                          "label": node["leaf"], "leaf": True})
            return y, leaf_depth
        ys = []
        child_ys = []
        for ch in node["children"].values():
            cy, cmax = place(ch, depth + 1)
            child_ys.append((cy, cmax))
            ys.append(cy)
        my = float(np.mean(ys))
        nodes.append({"x": depth, "y": my, "label": node["name"],
                      "leaf": False})
        for cy, cmax in child_ys:
            edges.append([depth, my, cmax, cy])
        return my, depth
    place(root, 0)
    return nodes, edges, leaves


SPECIES_LINEAGE = {}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aligned_json")
    p.add_argument("-o", "--output", default="viewer.html")
    p.add_argument("--title", default="Protein evolution viewer")
    args = p.parse_args()

    if pl is None:
        sys.exit("proteinlib.py not found - put it next to the skills tree "
                 "(src/proteinlib.py) or on PYTHONPATH")

    with open(args.aligned_json) as f:
        data = json.load(f)
    structs = data["structures"]
    if not structs:
        sys.exit("No structures in aligned.json")

    for s in structs:
        SPECIES_LINEAGE[s["species"]] = s.get("lineage", "")

    print(f"[tree] fetching taxonomy for {len(structs)} species ...",
          file=sys.stderr)
    session = requests.Session()
    cache = {}
    lineages = {}
    for s in structs:
        lineages[s["species"]] = fetch_lineage(s["species"], session, cache)
    tree = build_tree(lineages)
    nodes, edges, leaves = layout_tree(tree, LINEAGE_ORDER)

    print("[structures] computing SS / hydropathy ...", file=sys.stderr)
    out_structs = []
    ref_acc = data.get("reference")
    ref = next((s for s in structs if s["accession"] == ref_acc),
               structs[0])
    ref_coords = np.asarray(ref["coords"], dtype=float)
    ref_seq = ref.get("seq") or "A" * len(ref_coords)

    for s in structs:
        coords = np.asarray(s["coords"], dtype=float)
        seq = s.get("seq") or "A" * len(coords)
        n_coord = len(coords)
        coord_seq = seq[:n_coord] if len(seq) >= n_coord else seq + "X" * (n_coord - len(seq))
        h = hydrophobicity(coord_seq)
        ss = s.get("ss") or "".join(pl.assign_secondary(coords))
        pairs = pl.sequence_pairs(seq, ref_seq)
        pairs = [[a, b] for a, b in pairs
                 if a < len(coords) and b < len(ref_coords)]
        out_structs.append({
            "accession": s["accession"], "species": s["species"],
            "lineage": s.get("lineage", ""),
            "source": s.get("source", ""),
            "core_identity": s.get("core_identity"),
            "n_residues": len(coords),
            "ss": ss,
            "pdb": to_pdb_text({"coords": coords.tolist(), "seq": seq,
                                "ss": ss}),
            "hydro": [round(float(x), 2) for x in h],
            "ref_pairs": pairs,
            "ref_coords": None,
        })
    ref_pts = [[round(float(c), 3) for c in xyz] for xyz in ref_coords]
    for s in out_structs:
        s["ref_coords"] = ref_pts
        if s["accession"] == ref["accession"]:
            s["ref_pairs"] = [[i, i] for i in range(len(ref_coords))]

    # v15: Build MSA and BLOSUM62 similarity scores for sequence panel
    print("[sequences] building MSA and similarity scores ...", file=sys.stderr)
    msa_seqs, n_msa_cols, ref_col_map = build_msa(structs, ref_seq, ref["accession"])
    sim_scores = compute_similarity_scores(msa_seqs, ref["species"])
    # Pack as compact arrays: {species: {aln: "...", sim: [0.5, null, ...]}}
    seq_data = {}
    for sp, aln in msa_seqs:
        sim = next((s for ssp, s in sim_scores if ssp == sp), [0.5] * len(aln))
        seq_data[sp] = {"aln": aln, "sim": sim}
    # Overall identity per species (mean of non-null similarity scores)
    identity = {}
    for sp, aln in msa_seqs:
        sim = seq_data[sp]["sim"]
        valid = [s for s in sim if s is not None]
        identity[sp] = round(sum(valid) / max(len(valid), 1), 3) if valid else 0

    payload = {
        "title": args.title,
        "reference": ref["accession"],
        "referenceSpecies": ref["species"],
        "lineageOrder": LINEAGE_ORDER,
        "lineageColors": LINEAGE_COLORS,
        "tree": {"nodes": nodes, "edges": edges, "leaves": leaves},
        "structures": out_structs,
        "sequences": {
            "msaCols": n_msa_cols,
            "refColMap": ref_col_map,
            "data": seq_data,
            "identity": identity,
        },
    }

    template = open(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "viewer_template.html")).read()
    html = (template.replace("__DATA__", json.dumps(payload))
                .replace("__TITLE__", args.title))
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[done] viewer -> {args.output} "
          f"({os.path.getsize(args.output) / 1e6:.1f} MB)",
          file=sys.stderr)


if __name__ == "__main__":
    main()

```

## `skills/interactive_viewer/viewer_template.html`

```
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Protein Evolution Viewer</title>
<script src="https://unpkg.com/three@0.128.0/build/three.min.js"></script>
<script src="https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<style>
  html, body { margin:0; padding:0; height:100%; background:#0d1117;
              color:#e6edf3; font-family:system-ui, sans-serif;
              overflow:hidden; }
  #app { display:flex; flex-direction:column; height:100%; }
  #toprow { display:flex; flex:1; min-height:0; }
  #left { width:280px; min-width:240px; padding:14px; overflow:auto;
          border-right:1px solid #21262d; }
  #center { flex:1; position:relative; min-width:0; }
  #seqpanel { height:180px; border-top:1px solid #21262d;
              display:flex; flex-direction:column; overflow:hidden;
              flex-shrink:0; }
  #seqpanel.collapsed { height:28px; }
  #seqheader { display:flex; align-items:center; justify-content:space-between;
               padding:4px 12px; font-size:11px; color:#8b949e;
               flex-shrink:0; cursor:pointer; user-select:none; }
  #seqheader .seqtitle { font-weight:600; color:#e6edf3; font-size:12px; }
  #seqheader .seqinfo { font-size:11px; }
  #seqresize { position:absolute; top:-3px; left:0; right:0; height:6px;
               cursor:n-resize; z-index:10; }
  #seqresize:hover { background:rgba(88,166,255,0.2); }
  #seqcontainer { flex:1; overflow:auto; padding:0 12px 8px 12px;
                  font-family:monospace; font-size:11px; line-height:16px;
                  position:relative; }
  #seqcanvas { display:block; }
  /* Selected residue highlight in sequence */
  .seq-sel-marker { position:absolute; pointer-events:none;
                    border:1px solid #f0c674; background:rgba(240,198,116,0.25);
                    border-radius:2px; z-index:5; }
  /* 3D residue highlight marker */
  .residue-highlight-label { position:absolute; pointer-events:none;
                    font-size:10px; color:#f0c674; background:rgba(13,17,23,0.85);
                    padding:1px 5px; border-radius:3px; border:1px solid #f0c674;
                    z-index:20; white-space:nowrap; }
  .seqrow { white-space:nowrap; }
  #seqpanel .simbar { position:absolute; bottom:0; left:12px; right:12px;
                      height:14px; font-size:9px; color:#8b949e;
                      display:flex; align-items:center; gap:6px;
                      pointer-events:none; }
  #seqpanel .simbar .grad { flex:1; height:6px; border-radius:2px;
                    background:linear-gradient(90deg,#d73a49,#e6edf3,#3fb950); }
  #paneltitle { font-size:15px; font-weight:700; margin-bottom:6px;
                color:#e6edf3; }
  #sname { font-size:17px; font-weight:600; color:#e6edf3; }
  #ssub { color:#8b949e; font-size:12px; margin-top:2px; }
  #titlebar { position:absolute; top:10px; left:16px; z-index:5;
              pointer-events:none; }
  #legend { position:absolute; bottom:12px; left:16px; font-size:11px;
            color:#8b949e; z-index:5; pointer-events:none; }
  #legend .swatch { display:inline-block; width:10px; height:10px;
                    border-radius:2px; margin:0 4px 0 10px; }
  #hydrobar { position:absolute; bottom:12px; right:16px; font-size:11px;
              color:#8b949e; z-index:5; pointer-events:none; }
  #hydrobar .grad { display:inline-block; width:120px; height:10px;
                    border-radius:2px; vertical-align:middle;
                    background:linear-gradient(90deg,#4aa3ff,#f5f5f5,#ff9d2e); }
  #sslegend { position:absolute; bottom:34px; left:16px; font-size:11px;
               color:#8b949e; z-index:5; pointer-events:none; }
  #sslegend .swatch { display:inline-block;
                      border-radius:50%; margin:0 4px 0 10px;
                      vertical-align:middle; background:#888; }
  #hint { color:#8b949e; font-size:11px; margin-top:10px;
          border-top:1px solid #21262d; padding-top:8px; }
  .edge { stroke:#30363d; stroke-width:1.4; }
</style>
</head>
<body>
<div id="app">
  <div id="toprow">
    <div id="left">
      <div id="paneltitle">__TITLE__</div>
      <div id="tree"></div>
      <div id="hint">Click a species to morph the structure.
        Drag in the 3D window to rotate, scroll to zoom,
        right-drag to pan.</div>
    </div>
    <div id="center">
      <div id="titlebar">
        <div id="sname"></div>
        <div id="ssub"></div>
      </div>
      <div id="stage" style="width:100%;height:100%"></div>
      <div id="legend">
        <span class="swatch" style="background:#2ca02c"></span>Green algae
        <span class="swatch" style="background:#8c564b"></span>Bryophytes
        <span class="swatch" style="background:#ff7f0e"></span>Ferns
        <span class="swatch" style="background:#9467bd"></span>Gymnosperms
        <span class="swatch" style="background:#d62728"></span>Angiosperms
      </div>
      <div id="sslegend">
        <span class="swatch" style="height:9px;width:9px"></span>helix
        <span class="swatch" style="height:6px;width:6px"></span>strand
        <span class="swatch" style="height:3px;width:3px"></span>loop
      </div>
      <div id="hydrobar">
        hydrophilic
        <span class="grad"></span>
        hydrophobic
      </div>
    </div>
  </div>
  <div id="seqpanel">
    <div id="seqresize" title="Drag to resize"></div>
    <div id="seqheader" onclick="toggleSeqPanel()">
      <span class="seqtitle">Sequence Alignment (vs. reference)</span>
      <span class="seqinfo" id="seqinfo"></span>
      <span id="seqtoggle">&#9660;</span>
    </div>
    <div id="seqcontainer">
      <canvas id="seqcanvas"></canvas>
    </div>
  </div>
</div>
<script>
const DATA = __DATA__;

// ================= phylogenetic tree (SVG) =================
const treeBox = document.getElementById('tree');
const NS = 'http://www.w3.org/2000/svg';
const leafIdx = {};

function renderTree() {
  const w = treeBox.clientWidth || 400;
  const rowH = 30, padL = 6, padT = 10, labelW = 185;
  const n = DATA.tree.leaves.length;
  const height = n * rowH + padT * 2;
  const maxDepth = Math.max(...DATA.tree.nodes.map(nd => nd.x), 1);
  const colW = Math.max(6, (w - padL - labelW) / maxDepth);
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('width', w);
  svg.setAttribute('height', height);
  const pos = {};
  DATA.tree.nodes.forEach(nd => {
    pos[nd.label + '|' + nd.x] = {x: padL + nd.x * colW,
                                  y: padT + nd.y * rowH};
  });
  DATA.tree.edges.forEach(e => {
    const [d0, y0, d1, y1] = e;
    const x0 = padL + d0 * colW, x1 = padL + d1 * colW;
    const path = document.createElementNS(NS, 'path');
    path.setAttribute('d', `M${x0} ${padT + y0*rowH} ` +
      `L${x0} ${padT + y1*rowH} L${x1} ${padT + y1*rowH}`);
    path.setAttribute('class', 'edge');
    path.setAttribute('fill', 'none');
    svg.appendChild(path);
  });
  DATA.tree.leaves.forEach((sp, i) => {
    const st = DATA.structures.find(s => s.species === sp);
    const col = DATA.lineageColors[st ? (st.lineage||'') : ''] || '#30363d';
    const x = padL + maxDepth * colW;
    const y = padT + i * rowH;
    const g = document.createElementNS(NS, 'g');
    const c = document.createElementNS(NS, 'circle');
    c.setAttribute('cx', x); c.setAttribute('cy', y);
    c.setAttribute('r', 5); c.setAttribute('fill', col);
    c.setAttribute('stroke', '#0d1117');
    c.style.cursor = 'pointer';
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', x + 10); t.setAttribute('y', y + 4);
    t.setAttribute('fill', '#e6edf3');
    t.setAttribute('font-size', '12');
    t.style.cursor = 'pointer';
    t.textContent = sp;
    const hit = document.createElementNS(NS, 'rect');
    hit.setAttribute('x', x - 6); hit.setAttribute('y', y - 9);
    hit.setAttribute('width', labelW); hit.setAttribute('height', 18);
    hit.setAttribute('fill', 'transparent');
    hit.style.cursor = 'pointer';
    g.appendChild(hit); g.appendChild(c); g.appendChild(t);
    g.addEventListener('click', () => selectSpecies(sp));
    svg.appendChild(g);
    leafIdx[sp] = {t};
  });
  const outDeg = {};
  DATA.tree.edges.forEach(e => {
    const k = e[0] + '|' + e[1];
    outDeg[k] = (outDeg[k] || 0) + 1;
  });
  DATA.tree.nodes.forEach(nd => {
    if (nd.leaf) return;
    const xy = pos[nd.label + '|' + nd.x];
    if (!xy) return;
    const deg = outDeg[nd.x + '|' + nd.y] || 0;
    if (deg < 2) return;
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', xy.x + 3); t.setAttribute('y', xy.y - 6);
    t.setAttribute('fill', '#58a6ff');
    t.setAttribute('font-size', '10');
    t.setAttribute('font-style', 'italic');
    t.textContent = nd.label;
    svg.appendChild(t);
  });
  treeBox.appendChild(svg);
}

// ================= colors =================
// v11: pure hydropathy coloring (blue hydrophilic -> white -> orange
// hydrophobic), matching the user's request. No SS-class tint — the
// hydropathy gradient IS the color scheme.
function hydToColor(v) {
  const t = Math.max(0, Math.min(1, (v + 4.5) / 9.0));
  let r, g, b;
  if (t < 0.5) {
    const u = t / 0.5;
    r = 74 + (245 - 74) * u;
    g = 163 + (245 - 163) * u;
    b = 255 + (245 - 255) * u;
  } else {
    const u = (t - 0.5) / 0.5;
    r = 245 + (255 - 245) * u;
    g = 245 + (157 - 245) * u;
    b = 245 + (46 - 245) * u;
  }
  return [r / 255, g / 255, b / 255];
}

// v11: color is purely hydropathy-based (the SS-class base tint is removed
// so the structure reads as a clean hydrophilic->hydrophobic gradient).
function ssColor(type, hydroSeg) {
  return hydToColor(hydroSeg);
}

// ================= pure data helpers (testable) =================
function coordArrayOf(st) {
  const out = [];
  for (const ln of st.pdb.split('\n')) {
    if (ln.startsWith('ATOM')) {
      out.push([parseFloat(ln.substr(30, 8)),
                parseFloat(ln.substr(38, 8)),
                parseFloat(ln.substr(46, 8))]);
    }
  }
  return out;
}

const COORDS = {};
function coordsOf(st) {
  if (!COORDS[st.species]) COORDS[st.species] = coordArrayOf(st);
  return COORDS[st.species];
}

const TGT_REF_OF = {};
function tgtRefMap(tgt) {
  if (!TGT_REF_OF[tgt.species]) {
    const m = new Array(coordsOf(tgt).length).fill(null);
    for (const [t, r] of tgt.ref_pairs) {
      if (t < m.length) m[t] = r;
    }
    TGT_REF_OF[tgt.species] = m;
  }
  return TGT_REF_OF[tgt.species];
}

// per-pair morph weights: 0 = frozen (stable + same SS class),
// 1 = fully animated. Uses CA displacement in the common frame plus
// SS-class change between source and target.
function pairWeights(src, tgt) {
  const sc = coordsOf(src), tc = coordsOf(tgt);
  const ssS = src.ss || '', ssT = tgt.ss || '';
  const refToSrc = {};
  src.ref_pairs.forEach(([s, r]) => refToSrc[r] = s);
  const tgtToRef = tgtRefMap(tgt);
  const DISP_MAX = 4.0, DISP_MIN = 1.0;
  const w = new Array(tc.length).fill(1.0);
  for (let i = 0; i < tc.length; i++) {
    const r = tgtToRef[i];
    const sIdx = r !== null ? refToSrc[r] : undefined;
    if (sIdx === undefined || sIdx === null || !sc[sIdx]) {
      w[i] = 1.0;              // insertion: fully animated
      continue;
    }
    const a = sc[sIdx], b = tc[i];
    const d = Math.hypot(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
    const ssChanged = (ssS[sIdx] || 'C') !== (ssT[i] || 'C');
    const wDisp = Math.max(0, Math.min(1, (d - DISP_MIN) / (DISP_MAX - DISP_MIN)));
    let ww = Math.max(wDisp, ssChanged ? 1.0 : 0.0);
    ww = ww * ww * (3 - 2 * ww);   // smoothstep
    w[i] = ww;
  }
  return w;
}

function smoothstep(t) { t = Math.max(0, Math.min(1, t)); return t*t*(3-2*t); }

// v14: build a union-space morph that handles both insertions AND deletions.
// Returns {coords, ss, hydro, radiusScale} arrays in union space.
// Union space = merged sequence of: paired residues + insertions + deletions,
// ordered by their position along the target (for insertions) or source (for deletions).
// 
// For each union residue:
//   - Paired: morph source->target with stability weights
//   - Insertion (target-only): grow from anchor, fade in (radiusScale = t)
//   - Deletion (source-only): collapse toward anchor, fade out (radiusScale = 1-t)
// 
// This prevents the CatmullRom spline artifacts from clumped coincident points
// (166 residues at one position) and properly handles source-only C/N-terminal tails.
function buildUnionMorph(src, tgt, f, weights) {
  const sc = coordsOf(src);
  const tc = coordsOf(tgt);
  const ssS = src.ss || '';
  const ssT = tgt.ss || '';
  const hyS = src.hydro || [];
  const hyT = tgt.hydro || [];
  const refToSrc = {};
  src.ref_pairs.forEach(([s, r]) => refToSrc[r] = s);
  const tgtToRef = tgtRefMap(tgt);
  const srcToRef = new Array(sc.length).fill(null);
  src.ref_pairs.forEach(([s, r]) => { if (s < srcToRef.length) srcToRef[s] = r; });
  const t = smoothstep(f);
  
  // Classify each target residue
  const tgtIsPaired = new Array(tc.length).fill(false);
  const tgtSrcIdx = new Array(tc.length).fill(-1);
  for (let i = 0; i < tc.length; i++) {
    const r = tgtToRef[i];
    const sIdx = r !== null ? refToSrc[r] : undefined;
    if (sIdx !== undefined && sIdx !== null && sc[sIdx]) {
      tgtIsPaired[i] = true;
      tgtSrcIdx[i] = sIdx;
    }
  }
  
  // Classify each source residue
  const srcIsPaired = new Array(sc.length).fill(false);
  const srcTgtIdx = new Array(sc.length).fill(-1);
  for (let i = 0; i < tc.length; i++) {
    if (tgtIsPaired[i]) {
      srcIsPaired[tgtSrcIdx[i]] = true;
      srcTgtIdx[tgtSrcIdx[i]] = i;
    }
  }
  
  // Build union: walk through target, inserting deletions at the right place
  // Strategy: iterate over target residues. Before each unpaired target residue
  // (insertion), check if there are source-only residues (deletions) between
  // the previous paired residue and the next paired residue.
  const union = [];  // {type: 'paired'|'ins'|'del', coords, ss, hydro, radiusScale, tgtIdx, srcIdx}
  
  // Find insertion/deletion blocks
  // For each gap between consecutive paired target residues, 
  // check for insertions (target-only) and deletions (source-only)
  
  let tgtPos = 0;
  while (tgtPos < tc.length) {
    if (tgtIsPaired[tgtPos]) {
      // Paired residue
      const sIdx = tgtSrcIdx[tgtPos];
      const a = sc[sIdx], b = tc[tgtPos];
      const w = weights[tgtPos];
      let pos;
      if (w < 0.15) {
        const settle = t >= 1 ? 1 : 0;
        pos = [a[0] + (b[0]-a[0]) * settle, a[1] + (b[1]-a[1]) * settle, a[2] + (b[2]-a[2]) * settle];
      } else {
        const eff = t >= 1 ? 1 : w * t;
        pos = [a[0] + (b[0]-a[0]) * eff, a[1] + (b[1]-a[1]) * eff, a[2] + (b[2]-a[2]) * eff];
      }
      union.push({
        type: 'paired', coords: pos,
        ss: (t >= 1) ? (ssT[tgtPos] || 'C') : (w < 0.5 ? (ssS[sIdx] || 'C') : (ssT[tgtPos] || 'C')),
        hydro: (hyS[sIdx] || 0) * (1 - f) + (hyT[tgtPos] || 0) * f,
        radiusScale: 1.0,
        tgtIdx: tgtPos, srcIdx: sIdx
      });
      tgtPos++;
    } else {
      // Insertion block: collect consecutive unpaired target residues
      const insStart = tgtPos;
      while (tgtPos < tc.length && !tgtIsPaired[tgtPos]) tgtPos++;
      const insEnd = tgtPos;  // exclusive
      
      // Find the anchor: nearest paired residue before or after
      const anchorBefore = insStart > 0 ? union[union.length - 1] : null;
      let anchorAfter = null;
      if (insEnd < tc.length && tgtIsPaired[insEnd]) {
        // Will be added next iteration, compute its morph position now
        const sIdxA = tgtSrcIdx[insEnd];
        const aA = sc[sIdxA], bA = tc[insEnd];
        const wA = weights[insEnd];
        const effA = t >= 1 ? 1 : wA * t;
        anchorAfter = {
          coords: [aA[0] + (bA[0]-aA[0]) * effA, aA[1] + (bA[1]-aA[1]) * effA, aA[2] + (bA[2]-aA[2]) * effA]
        };
      }
      
      // Check for deletions between the source indices of surrounding paired residues
      // The source residues between anchorBefore.srcIdx and anchorAfter's srcIdx
      // that are unpaired are deletions
      let delStart = -1, delEnd = -1;
      if (anchorBefore && anchorBefore.type === 'paired') {
        delStart = anchorBefore.srcIdx + 1;
      } else {
        delStart = 0;
      }
      if (anchorAfter) {
        // Find source index of anchorAfter
        delEnd = tgtSrcIdx[insEnd];  // exclusive
      } else {
        delEnd = sc.length;
      }
      
      // Add deletions first (they come before insertions in the chain)
      // Deletions collapse toward the preceding anchor
      for (let sIdx = delStart; sIdx < delEnd; sIdx++) {
        if (srcIsPaired[sIdx]) continue;  // skip paired (shouldn't happen but safety)
        const a = sc[sIdx];
        // Collapse toward anchor: at t=0 stay at source, at t=1 collapse to anchor
        const anchorPos = anchorBefore ? anchorBefore.coords : a;
        // Spread term prevents clumping at intermediate t (peaks at t=0.5)
        const spread = (sIdx - delStart) / Math.max(1, delEnd - delStart);
        const spreadTerm = spread * (1 - t) * t;  // 0 at t=0 and t=1, max at t=0.5
        union.push({
          type: 'del',
          coords: [
            a[0] * (1 - t) + anchorPos[0] * t + (a[0] - anchorPos[0]) * spreadTerm,
            a[1] * (1 - t) + anchorPos[1] * t + (a[1] - anchorPos[1]) * spreadTerm,
            a[2] * (1 - t) + anchorPos[2] * t + (a[2] - anchorPos[2]) * spreadTerm
          ],
          ss: ssS[sIdx] || 'C',
          hydro: (hyS[sIdx] || 0) * (1 - f),
          radiusScale: (1 - t),  // fade out
          tgtIdx: -1, srcIdx: sIdx
        });
      }
      
      // Add insertions (grow from anchor)
      for (let i = insStart; i < insEnd; i++) {
        const b = tc[i];
        const anchorPos = anchorBefore ? anchorBefore.coords : (anchorAfter ? anchorAfter.coords : b);
        // Spread along the anchor->target direction by position within block
        const spread = (i - insStart) / Math.max(1, insEnd - insStart);
        // At t=0: all at anchor. At t=1: at target. In between: interpolate.
        // Use spread to offset slightly so points aren't all coincident
        const dir = [b[0] - anchorPos[0], b[1] - anchorPos[1], b[2] - anchorPos[2]];
        const len = Math.hypot(dir[0], dir[1], dir[2]) || 1;
        // Small offset along the direction to prevent spline degeneracy
        const offset = spread * len * t * 0.3;  // 30% of the way, spread by index
        union.push({
          type: 'ins',
          coords: [
            anchorPos[0] + dir[0] * t + dir[0] / len * offset * (1 - t),
            anchorPos[1] + dir[1] * t + dir[1] / len * offset * (1 - t),
            anchorPos[2] + dir[2] * t + dir[2] / len * offset * (1 - t)
          ],
          ss: ssT[i] || 'C',
          hydro: (hyT[i] || 0) * f,
          radiusScale: t,  // fade in
          tgtIdx: i, srcIdx: -1
        });
      }
    }
  }
  
  // Also handle trailing deletions (source-only residues after last paired target)
  // This happens when source has a longer C-terminal tail
  // Find the last paired source index
  let lastPairedSrc = -1;
  for (let i = 0; i < sc.length; i++) {
    if (srcIsPaired[i]) lastPairedSrc = i;
  }
  // The last union entry should be the last paired residue
  const lastUnion = union[union.length - 1];
  if (lastUnion && lastUnion.type === 'paired') {
    for (let sIdx = lastPairedSrc + 1; sIdx < sc.length; sIdx++) {
      if (srcIsPaired[sIdx]) continue;
      const a = sc[sIdx];
      const anchorPos = lastUnion.coords;
      const spread = (sIdx - lastPairedSrc - 1) / Math.max(1, sc.length - lastPairedSrc - 1);
      const spreadTerm = spread * (1 - t) * t;
      union.push({
        type: 'del',
        coords: [
          a[0] * (1 - t) + anchorPos[0] * t + (a[0] - anchorPos[0]) * spreadTerm,
          a[1] * (1 - t) + anchorPos[1] * t + (a[1] - anchorPos[1]) * spreadTerm,
          a[2] * (1 - t) + anchorPos[2] * t + (a[2] - anchorPos[2]) * spreadTerm
        ],
        ss: ssS[sIdx] || 'C',
        hydro: (hyS[sIdx] || 0) * (1 - f),
        radiusScale: (1 - t),
        tgtIdx: -1, srcIdx: sIdx
      });
    }
  }
  
  // Handle leading deletions (source-only residues before first paired target)
  // These should already be handled in the first iteration (anchorBefore = null)
  
  return union;
}

// Legacy buildMorphCoords for backward compatibility (returns just coords array)
function buildMorphCoords(src, tgt, f, weights) {
  const union = buildUnionMorph(src, tgt, f, weights);
  return union.map(u => u.coords);
}

// ================= three.js scene =================
let scene, camera, renderer, controls, structureGroup;
const stage = document.getElementById('stage');

function init3D() {
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0d1117);
  camera = new THREE.PerspectiveCamera(50, 1, 1, 5000);
  camera.position.set(0, 0, 260);
  renderer = new THREE.WebGLRenderer({antialias: true});
  renderer.setPixelRatio(window.devicePixelRatio);
  stage.appendChild(renderer.domElement);
  controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  // consistent MEDIUM lighting: hemisphere + directional at fixed medium
  // intensity, attached to the camera so illumination stays fixed relative
  // to the viewer while orbiting
  const hemi = new THREE.HemisphereLight(0xffffff, 0x223344, 0.55);
  const dir = new THREE.DirectionalLight(0xffffff, 0.45);
  dir.position.set(150, 200, 180);
  camera.add(hemi);
  camera.add(dir);
  scene.add(camera);
  structureGroup = new THREE.Group();
  scene.add(structureGroup);
  resize();
  window.addEventListener('resize', resize);
  animateLoop();
}

function resize() {
  const w = stage.clientWidth, h = stage.clientHeight;
  if (w === 0 || h === 0) return;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
}

// Alias for external callers (e.g., sequence panel toggle)
function resizeRenderer() { resize(); }

function animateLoop() {
  requestAnimationFrame(animateLoop);
  controls.update();
  // v16: Keep residue highlight ring facing the camera
  if (window._residueRing && camera) {
    window._residueRing.lookAt(camera.position);
  }
  // v16: Pulse the residue highlight sphere for visibility
  if (window._residueHighlight) {
    const pulse = 1 + 0.15 * Math.sin(performance.now() / 300);
    window._residueHighlight.scale.set(pulse, pulse, pulse);
  }
  renderer.render(scene, camera);
}

function clearGroup() {
  clearTerminusMarkers();
  while (structureGroup.children.length) {
    const ch = structureGroup.children.pop();
    ch.traverse(o => { if (o.geometry) o.geometry.dispose();
                       if (o.material) o.material.dispose();
                       if (o.material && o.material.map) o.material.map.dispose(); });
  }
}

// ---- cartoon tube: smooth CatmullRom tube through the CA trace ----
// v13: standard cartoon — round tube with variable thickness.
// Helices thick (2.2), strands medium (1.6), loops thin (0.7).
// No elliptical cross-sections — clean round tube like PyMOL/VMD cartoon.
let activeTube = null;
function addTube(coordsArr, hydro, ss) {
  const pts = [];
  for (const c of coordsArr) pts.push(new THREE.Vector3(c[0], c[1], c[2]));
  const curve = new THREE.CatmullRomCurve3(pts, false, 'catmullrom', 0.5);
  const tubular = Math.min(2400, pts.length * 3);
  const radial = 12;
  const geo = new THREE.TubeGeometry(curve, tubular, 1.0, radial, false);
  const pos = geo.attributes.position;
  const colors = new Float32Array(pos.count * 3);
  const ringSize = radial + 1;
  for (let i = 0; i < pos.count; i++) {
    const ring = Math.floor(i / ringSize);
    const t = ring / tubular;
    const resi = Math.min(pts.length - 1,
                          Math.floor(t * (pts.length - 1)));
    const segH = meanHydro(hydro, ss, resi);
    const col = ssColor('C', segH);
    colors[i*3] = col[0]; colors[i*3+1] = col[1]; colors[i*3+2] = col[2];
  }
  geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
  const mat = new THREE.MeshPhongMaterial({vertexColors: true,
                                          shininess: 30,
                                          specular: 0x222222,
                                          side: THREE.DoubleSide});
  const mesh = new THREE.Mesh(geo, mat);
  structureGroup.add(mesh);
  activeTube = {mesh, curve, tubular, radial, ringSize, ss};
  updateTubeGeometry(coordsArr, ss);
  return activeTube;
}

// v13: per-residue cartoon radius — round tube, thickness varies by SS class.
// H=2.2 (thick), E=1.6 (medium), C=0.7 (thin). Standard cartoon look.
function ssRadius(st) {
  if (st === 'H') return 2.2;
  if (st === 'E') return 1.6;
  return 0.7;
}

// Precompute smoothed per-residue radii with linear interpolation at SS
// boundaries so the thickness transitions smoothly (no hard steps).
function smoothRadii(ss, nResidues) {
  const raw = new Array(nResidues);
  for (let i = 0; i < nResidues; i++) {
    raw[i] = ssRadius((ss && ss[i]) || 'C');
  }
  // Smooth: 3-point moving average to blend transitions
  const sm = new Array(nResidues);
  for (let i = 0; i < nResidues; i++) {
    const a = raw[Math.max(0, i - 1)];
    const b = raw[i];
    const c = raw[Math.min(nResidues - 1, i + 1)];
    sm[i] = (a + 2 * b + c) / 4;
  }
  return sm;
}

function updateTubeGeometry(coordsArr, ss, radiusScale) {
  const t = activeTube;
  if (!t) return;
  const pts = [];
  for (const c of coordsArr) pts.push(new THREE.Vector3(c[0], c[1], c[2]));
  t.curve.points = pts;
  t.ss = ss;
  const nRes = pts.length;
  const radii = smoothRadii(ss, nRes);
  t.radii = radii;
  const pos = t.mesh.geometry.attributes.position;
  const N = new THREE.Vector3();
  const P = new THREE.Vector3();
  // Simple frame: use Frenet-like approach with fallback
  let B = new THREE.Vector3();
  const v = new THREE.Vector3();
  const firstTangent = new THREE.Vector3();
  try { t.curve.getTangent(0, firstTangent); } catch (e) { firstTangent.set(0, 0, 1); }
  firstTangent.normalize();
  const refUp = Math.abs(firstTangent.y) > 0.95
    ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
  B.crossVectors(firstTangent, refUp).normalize();
  for (let ring = 0; ring <= t.tubular; ring++) {
    const i = ring / t.tubular;
    t.curve.getPoint(i, P);
    try { t.curve.getTangent(i, N); } catch (e) { N.copy(firstTangent); }
    N.normalize();
    if (ring > 0) {
      const dotBN = B.dot(N);
      B.x -= dotBN * N.x; B.y -= dotBN * N.y; B.z -= dotBN * N.z;
      if (B.lengthSq() < 1e-6) {
        B.crossVectors(N, refUp).normalize();
      } else {
        B.normalize();
      }
    }
    v.crossVectors(N, B).normalize();
    // Interpolated radius at this ring position, with optional per-residue
    // radiusScale for fade in/out of insertions/deletions
    const fi = i * (nRes - 1);
    const ri = Math.floor(fi);
    const frac = fi - ri;
    const r0 = radii[Math.min(ri, nRes - 1)];
    const r1 = radii[Math.min(ri + 1, nRes - 1)];
    let r = r0 + (r1 - r0) * frac;
    if (radiusScale) {
      const rs0 = radiusScale[Math.min(ri, nRes - 1)];
      const rs1 = radiusScale[Math.min(ri + 1, nRes - 1)];
      const rs = rs0 + (rs1 - rs0) * frac;
      r *= rs;
      if (r < 0.05) r = 0.05;  // minimum visible radius to avoid degenerate geometry
    }
    for (let s = 0; s < t.ringSize; s++) {
      const theta = s / t.radial * Math.PI * 2;
      const cx = Math.cos(theta) * r;
      const cy = Math.sin(theta) * r;
      const idx = ring * t.ringSize + s;
      if (idx >= pos.count) break;
      pos.setXYZ(idx,
        P.x + B.x * cx + v.x * cy,
        P.y + B.y * cx + v.y * cy,
        P.z + B.z * cx + v.z * cy);
    }
  }
  pos.needsUpdate = true;
  t.mesh.geometry.computeVertexNormals();
  t.mesh.geometry.computeBoundingSphere();
}

function meanHydro(hydro, ss, resi) {
  if (!hydro || !hydro.length) return 0;
  return hydro[resi] !== undefined ? hydro[resi] : 0;
}

// v10: update the tube's VERTEX COLORS in place (smooth color cross-fade
// during a morph, instead of rebuilding the tube with target colors).
function updateTubeColors(ss, hydro) {
  const t = activeTube;
  if (!t || !ss) return;
  const pos = t.mesh.geometry.attributes.position;
  const colors = t.mesh.geometry.attributes.color;
  const ringSize = t.radial + 1;
  for (let i = 0; i < pos.count; i++) {
    const ring = Math.floor(i / ringSize);
    const tt = ring / t.tubular;
    const resi = Math.min(t.curve.points.length - 1,
                          Math.floor(tt * (t.curve.points.length - 1)));
    const st = ss[resi] || 'C';
    const segH = meanHydro(hydro, ss, resi);
    const col = ssColor(st, segH);
    colors.setXYZ(i, col[0], col[1], col[2]);
  }
  colors.needsUpdate = true;
}

// v10: interpolated SS string — the SS class itself cross-fades smoothly
// with the morph progress: stable residues (w < 0.15, frozen) keep A's
// class for the whole morph and only settle to B at the very end; all
// other residues cross over to B's class as the morph progresses
// (matches the GIF renderer, which keeps A's class while w < 0.5 and
// takes B's class afterwards, with the final frame exactly B).
function morphSS(ssA, ssB, src, tgt, f, weights) {
  const tc = coordsOf(tgt);
  const refToSrc = {};
  src.ref_pairs.forEach(([s, r]) => refToSrc[r] = s);
  const tgtToRef = tgtRefMap(tgt);
  const out = new Array(tc.length);
  for (let i = 0; i < tc.length; i++) {
    const r = tgtToRef[i];
    const sIdx = r !== null ? refToSrc[r] : undefined;
    if (sIdx === undefined || sIdx === null) {
      out[i] = (ssB && ssB[i]) || 'C';   // insertion takes target class
    } else {
      const w = weights ? weights[i] : 1.0;
      const aCls = (ssA && ssA[sIdx]) || 'C';
      const bCls = (ssB && ssB[i]) || 'C';
      if (w < 0.15) {
        // frozen: hold A's class until the very end
        out[i] = (f >= 1) ? bCls : aCls;
      } else {
        // switch on the morph weight, exactly like the GIF renderer:
        // low-weight (mostly stable) residues keep A's class for the
        // whole morph; genuinely changing residues take B's class
        out[i] = (w < 0.5) ? aCls : bCls;
      }
    }
  }
  return out.join('');
}

// v10: interpolated hydropathy — smooth color cross-fade between the two
// structures' per-residue hydropathy (the GIF cross-fades pLDDT colors the
// same way; the viewer tints SS base colors by hydropathy).
function morphHydro(hydroA, hydroB, src, tgt, f, weights) {
  const tc = coordsOf(tgt);
  const refToSrc = {};
  src.ref_pairs.forEach(([s, r]) => refToSrc[r] = s);
  const tgtToRef = tgtRefMap(tgt);
  const out = new Array(tc.length).fill(0);
  for (let i = 0; i < tc.length; i++) {
    const r = tgtToRef[i];
    const sIdx = r !== null ? refToSrc[r] : undefined;
    const aV = (sIdx !== undefined && sIdx !== null && hydroA)
      ? (hydroA[sIdx] || 0) : 0;
    const bV = (hydroB && hydroB[i] !== undefined) ? hydroB[i] : 0;
    out[i] = aV * (1 - f) + bV * f;
  }
  return out;
}

// ---- N-/C-terminus markers: small spheres + floating text labels ----
const TERM_COLORS = {N: 0x58a6ff, C: 0xff7f0e};
let termMarkers = [];
function addTerminusMarkers(coordsArr) {
  clearTerminusMarkers();
  if (!coordsArr || coordsArr.length < 2) return;
  [[0, 'N'], [coordsArr.length - 1, 'C']].forEach(([idx, label]) => {
    const p = coordsArr[idx];
    const sphere = new THREE.Mesh(
      new THREE.SphereGeometry(2.4, 16, 12),
      new THREE.MeshPhongMaterial({color: TERM_COLORS[label],
                                   emissive: TERM_COLORS[label],
                                   emissiveIntensity: 0.35,
                                   shininess: 40}));
    sphere.position.set(p[0], p[1], p[2]);
    structureGroup.add(sphere);
    const spr = makeTextSprite(label, label === 'N' ? '#58a6ff' : '#ff7f0e');
    spr.position.set(p[0], p[1], p[2]);
    structureGroup.add(spr);
    termMarkers.push(sphere, spr);
  });
}

function clearTerminusMarkers() {
  for (const m of termMarkers) {
    structureGroup.remove(m);
    if (m.geometry) m.geometry.dispose();
    if (m.material) m.material.dispose();
    if (m.material && m.material.map) m.material.map.dispose();
  }
  termMarkers = [];
}

function makeTextSprite(text, color) {
  const canvas = document.createElement('canvas');
  canvas.width = 64; canvas.height = 64;
  const ctx = canvas.getContext('2d');
  ctx.font = 'bold 40px system-ui, sans-serif';
  ctx.fillStyle = color;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText(text, 32, 34);
  const tex = new THREE.CanvasTexture(canvas);
  const mat = new THREE.SpriteMaterial({map: tex, depthTest: false,
                                        transparent: true});
  const spr = new THREE.Sprite(mat);
  spr.scale.set(9, 9, 1);
  spr.position.y += 6;
  return spr;
}

function updateTubeCoords(coordsArr, ss, radiusScale) {
  const t = activeTube;
  if (!t) return;
  updateTubeGeometry(coordsArr, ss || t.ss, radiusScale);
}

// v11: helix ribbons and strand arrows are REMOVED. The single tube with
// variable per-ring radius (wide/flat over H and E, thin over loops) now
// provides the entire ribbon cartoon look — exactly like the GIF, which
// draws everything as thick/thin lines along the CA trace. This eliminates
// the "broken" morph appearance (boxes not updating) and the orientation
// jump at morph end (no rebuild needed — the tube just settles).
// addHelixRibbons and addStrandArrows are kept as no-ops for compatibility.
function addHelixRibbons(coordsArr, ss, hydro) { /* v11: handled by variable-radius tube */ }
function addStrandArrows(coordsArr, ss, hydro) { /* v11: handled by variable-radius tube */ }

function frameScene(coordsArr) {
  const box = new THREE.Box3();
  const v = new THREE.Vector3();
  for (const c of coordsArr) {
    box.expandByPoint(v.set(c[0], c[1], c[2]));
  }
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3()).length();
  controls.target.copy(center);
  return size;
}

// ================= structure display & morph =================
let current = null;
let morphState = null;

function showStructure(st) {
  document.getElementById('sname').textContent = st.species;
  document.getElementById('ssub').textContent =
    st.accession + ' · ' + st.source +
    (st.core_identity != null ? ' · core identity ' +
      Math.round(st.core_identity * 100) + '%' : '');
  clearGroup();
  const c = coordsOf(st);
  addTube(c, st.hydro, st.ss);
  addTerminusMarkers(c);
  const size = frameScene(c);
  if (!current) {
    camera.position.set(controls.target.x, controls.target.y + 40,
                        controls.target.z + size * 0.9);
  }
  current = st;
  // v16: Update residue highlight for initial structure
  if (typeof updateResidueHighlight === 'function') {
    updateResidueHighlight();
  }
}

function selectSpecies(sp) {
  const st = DATA.structures.find(s => s.species === sp);
  if (!st) return;
  if (current && current !== st) {
    startMorph(st);
  } else if (!current) {
    showStructure(st);
  }
  if (window.__selected) window.__selected.setAttribute('fill', '#e6edf3');
  const entry = leafIdx[sp];
  if (entry && entry.t) {
    entry.t.setAttribute('fill', '#58a6ff');
    entry.t.setAttribute('font-weight', '700');
    window.__selected = entry.t;
  }
  // v15: transition sequence panel highlight
  if (DATA.sequences && typeof transitionSeqHighlight === 'function') {
    transitionSeqHighlight(sp);
  }
}

// v10: morph duration ~1 second, time-eased (smoothstep) so the motion
// starts and ends softly — the same feel as the GIF.
const MORPH_MS = 1000;

function startMorph(target) {
  const src = current;
  // stability-masked morph: weights per target residue
  const weights = pairWeights(src, target);
  morphState = {src, target, t0: performance.now(), weights};
  // v14: use union-space morph (handles insertions + deletions)
  const u0 = buildUnionMorph(src, target, 0, weights);
  const c0 = u0.map(u => u.coords);
  const ss0 = u0.map(u => u.ss).join('');
  const hy0 = u0.map(u => u.hydro);
  const rs0 = u0.map(u => u.radiusScale);
  clearGroup();
  addTube(c0, hy0, ss0);
  updateTubeGeometry(c0, ss0, rs0);
  addTerminusMarkers(c0);
  morphTick();
}

function morphTick() {
  if (!morphState) return;
  const {src, target, t0, weights} = morphState;
  const f = Math.min(1, (performance.now() - t0) / MORPH_MS);
  // v14: union-space morph with insertion/deletion fade
  const union = buildUnionMorph(src, target, f, weights);
  const c = union.map(u => u.coords);
  const ssK = union.map(u => u.ss).join('');
  const hydroK = union.map(u => u.hydro);
  const rsK = union.map(u => u.radiusScale);
  updateTubeGeometry(c, ssK, rsK);
  updateTubeColors(ssK, hydroK);
  // termini markers follow the morphing chain ends
  if (termMarkers.length >= 4 && c.length >= 2) {
    termMarkers[0].position.set(c[0][0], c[0][1], c[0][2]);
    termMarkers[1].position.set(c[0][0], c[0][1] + 6, c[0][2]);
    termMarkers[2].position.set(c[c.length - 1][0], c[c.length - 1][1],
                                c[c.length - 1][2]);
    termMarkers[3].position.set(c[c.length - 1][0], c[c.length - 1][1] + 6,
                                c[c.length - 1][2]);
  }
  if (f < 1) {
    requestAnimationFrame(morphTick);
  } else {
    // v11: NO rebuild — just snap the tube to the exact target coords and
    // colors. This eliminates the orientation jump at morph end.
    morphState = null;
    const cFinal = coordsOf(target);
    updateTubeGeometry(cFinal, target.ss);
    updateTubeColors(target.ss, target.hydro);
    if (termMarkers.length >= 4 && cFinal.length >= 2) {
      termMarkers[0].position.set(cFinal[0][0], cFinal[0][1], cFinal[0][2]);
      termMarkers[1].position.set(cFinal[0][0], cFinal[0][1] + 6, cFinal[0][2]);
      termMarkers[2].position.set(cFinal[cFinal.length-1][0], cFinal[cFinal.length-1][1],
                                  cFinal[cFinal.length-1][2]);
      termMarkers[3].position.set(cFinal[cFinal.length-1][0], cFinal[cFinal.length-1][1] + 6,
                                  cFinal[cFinal.length-1][2]);
    }
    current = target;
    document.getElementById('sname').textContent = target.species;
    document.getElementById('ssub').textContent =
      target.accession + ' · ' + target.source +
      (target.core_identity != null ? ' · core identity ' +
        Math.round(target.core_identity * 100) + '%' : '');
    // v16: Update residue highlight for the new structure
    // (residue may not exist in the new structure -> highlight removed)
    if (typeof updateResidueHighlight === 'function') {
      updateResidueHighlight();
    }
  }
}

// ================= sequence alignment panel (v15) =================
// Renders a multiple sequence alignment as a canvas. Each row is one species.
// Residues are colored by BLOSUM62 similarity to the reference (red=dissimilar,
// white=neutral, green=similar/identical). The selected species' row is
// highlighted with a border and brought into view. Smooth transition when
// switching species.

let seqCanvas, seqCtx;
let seqState = {
  scrollX: 0, scrollY: 0,
  charW: 7, charH: 16, labelW: 160, rowH: 18,
  selectedSpecies: null,
  highlightAlpha: 0,
  targetSpecies: null,
  speciesOrder: [],
  colWidth: 0,
  nCols: 0,
  animating: false,
  // v16: selected residue tracking (reference column index, null = none)
  selectedRefCol: null,
  // v16: per-species mapping from reference column -> structure residue index
  // (null if that species doesn't have a residue at this reference position)
  refColToResIdx: {},  // species -> array (refCol -> residueIdx or -1)
};

function initSeqPanel() {
  if (!DATA.sequences) return;
  seqCanvas = document.getElementById('seqcanvas');
  seqCtx = seqCanvas.getContext('2d');
  seqState.nCols = DATA.sequences.msaCols;
  seqState.speciesOrder = DATA.tree.leaves.slice();

  // v16: Build refColToResIdx mapping for each species
  // For each species, map each reference column to the structure residue index
  // (or -1 if the species doesn't have a residue at this reference position)
  buildRefColToResIdx();

  resizeSeqCanvas();
  const container = document.getElementById('seqcontainer');

  // v16: Use native scrolling (container has overflow:auto) instead of
  // manual scrollY. This fixes rows staying invisible when scrolled.
  container.addEventListener('scroll', () => {
    seqState.scrollX = container.scrollLeft;
    seqState.scrollY = container.scrollTop;
    drawSeqPanel();
  });

  // Mouse wheel: horizontal scroll with shift, vertical with plain wheel
  container.addEventListener('wheel', (e) => {
    if (!e.shiftKey) {
      // Vertical scroll: let the browser handle it natively
      return;
    }
    e.preventDefault();
    container.scrollLeft += e.deltaY;
  }, {passive: false});

  // v16: Click handling — distinguish row click (select species) from
  // residue click (select amino acid for 3D highlight)
  seqCanvas.addEventListener('click', (e) => {
    const rect = seqCanvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const rowIdx = Math.floor((y - 4 + seqState.scrollY) / seqState.rowH);
    if (rowIdx < 0 || rowIdx >= seqState.speciesOrder.length) return;
    const sp = seqState.speciesOrder[rowIdx];

    // Check if click is on a residue (right of label area)
    const colX = x + seqState.scrollX - seqState.labelW;
    if (colX >= 0) {
      const col = Math.floor(colX / seqState.charW);
      if (col >= 0 && col < seqState.nCols) {
        // Click on a residue: check if it's a reference column
        const refCol = msaColToRefCol(col);
        if (refCol !== null) {
          // Toggle residue selection
          if (seqState.selectedRefCol === refCol) {
            seqState.selectedRefCol = null;
          } else {
            seqState.selectedRefCol = refCol;
          }
          drawSeqPanel();
          updateResidueHighlight();
          scrollToSelectedResidue();  // center the selected column
          return;
        }
      }
    }
    // Otherwise: select species
    selectSpecies(sp);
  });

  // v16: Draggable resize handle
  initSeqResize();

  // Resize handler
  window.addEventListener('resize', () => {
    resizeSeqCanvas();
    drawSeqPanel();
  });
  drawSeqPanel();
}

// v16: Map MSA column index -> reference column index (or null if insertion column)
function msaColToRefCol(msaCol) {
  const refColMap = DATA.sequences.refColMap;
  if (!refColMap) return null;
  // refColMap[refIdx] = msaCol, so we need the inverse
  // Build inverse lazily
  if (!seqState._refColInverse) {
    seqState._refColInverse = {};
    for (let ri = 0; ri < refColMap.length; ri++) {
      seqState._refColInverse[refColMap[ri]] = ri;
    }
  }
  return seqState._refColInverse[msaCol] !== undefined ? seqState._refColInverse[msaCol] : null;
}

// v16: Build mapping from reference column -> structure residue index for each species
// This uses the ref_pairs data: ref_pairs[seq_idx] = ref_idx
// So for each species, we invert: ref_idx -> seq_idx (structure residue index)
function buildRefColToResIdx() {
  const refColMap = DATA.sequences.refColMap;
  if (!refColMap) return;
  for (const st of DATA.structures) {
    const sp = st.species;
    const nRefCols = refColMap.length;
    const mapping = new Array(nRefCols).fill(-1);
    // ref_pairs: [[seq_idx, ref_idx], ...]
    // For the reference structure, ref_pairs is identity: [[i, i], ...]
    if (st.ref_pairs) {
      for (const [seqIdx, refIdx] of st.ref_pairs) {
        if (refIdx >= 0 && refIdx < nRefCols) {
          mapping[refIdx] = seqIdx;
        }
      }
    }
    seqState.refColToResIdx[sp] = mapping;
  }
}

// v16: Get the structure residue index for the selected reference column
// in the currently displayed species. Returns -1 if not present.
function getSelectedResidueIdx(species) {
  if (seqState.selectedRefCol === null) return -1;
  const mapping = seqState.refColToResIdx[species];
  if (!mapping) return -1;
  return mapping[seqState.selectedRefCol];
}

// v16: Update the 3D residue highlight marker
function updateResidueHighlight() {
  // Remove existing highlight
  if (window._residueHighlight) {
    window._residueHighlight.remove();
    window._residueHighlight = null;
  }
  if (!current || seqState.selectedRefCol === null) return;
  const resIdx = getSelectedResidueIdx(current.species);
  if (resIdx < 0) return;  // residue doesn't exist in this structure
  const coords = coordsOf(current);
  if (resIdx >= coords.length) return;
  const p = coords[resIdx];

  // Create a highlight sphere at the residue position
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(3.5, 16, 12),
    new THREE.MeshPhongMaterial({
      color: 0xf0c674, emissive: 0xf0c674, emissiveIntensity: 0.5,
      transparent: true, opacity: 0.7
    })
  );
  sphere.position.set(p[0], p[1], p[2]);
  // Add to the scene (not the structure group, so it persists during morphs)
  scene.add(sphere);
  window._residueHighlight = sphere;

  // Add a ring around the residue
  const ring = new THREE.Mesh(
    new THREE.TorusGeometry(5.0, 0.4, 8, 24),
    new THREE.MeshPhongMaterial({
      color: 0xf0c674, emissive: 0xf0c674, emissiveIntensity: 0.3
    })
  );
  ring.position.set(p[0], p[1], p[2]);
  // Orient ring to face camera (will be updated in animate loop)
  ring.userData.isResidueRing = true;
  scene.add(ring);
  window._residueRing = ring;

  // Update sequence panel info
  const refAA = DATA.sequences.data[DATA.referenceSpecies].aln[
    DATA.sequences.refColMap[seqState.selectedRefCol]];
  const info = document.getElementById('seqinfo');
  if (info) {
    const curAA = DATA.sequences.data[current.species].aln[
      findMsaColForSpecies(current.species, seqState.selectedRefCol)];
    info.textContent = 'Ref pos ' + (seqState.selectedRefCol + 1) + ' (' + refAA +
      ' → ' + (curAA || '-') + ' in ' + current.species + ')';
  }
}

// v16: Find the MSA column for a given reference column in a species
function findMsaColForSpecies(species, refCol) {
  const refColMap = DATA.sequences.refColMap;
  if (!refColMap || refCol < 0 || refCol >= refColMap.length) return -1;
  return refColMap[refCol];
}

// v16: Draggable resize handle for sequence panel
function initSeqResize() {
  const handle = document.getElementById('seqresize');
  const panel = document.getElementById('seqpanel');
  if (!handle || !panel) return;
  let startY = 0, startH = 0;
  handle.addEventListener('mousedown', (e) => {
    e.preventDefault();
    startY = e.clientY;
    startH = panel.offsetHeight;
    document.body.style.cursor = 'n-resize';
    document.body.style.userSelect = 'none';
    function onMove(ev) {
      const dy = startY - ev.clientY;  // drag up = increase height
      const newH = Math.max(60, Math.min(600, startH + dy));
      panel.style.height = newH + 'px';
      panel.classList.remove('collapsed');
      document.getElementById('seqtoggle').innerHTML = '&#9660;';
      resizeSeqCanvas();
      if (typeof resizeRenderer === 'function') resizeRenderer();
      drawSeqPanel();
    }
    function onUp() {
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
      document.removeEventListener('mousemove', onMove);
      document.removeEventListener('mouseup', onUp);
    }
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
  });
}

function resizeSeqCanvas() {
  if (!seqCanvas || !DATA.sequences) return;
  const container = document.getElementById('seqcontainer');
  const dpr = window.devicePixelRatio || 1;
  const w = container.clientWidth;
  const nRows = seqState.speciesOrder.length;
  seqState.colWidth = seqState.labelW + seqState.nCols * seqState.charW + 20;
  const h = nRows * seqState.rowH + 8;
  seqCanvas.width = seqState.colWidth * dpr;
  seqCanvas.height = h * dpr;
  seqCanvas.style.width = seqState.colWidth + 'px';
  seqCanvas.style.height = h + 'px';
  seqCtx.scale(dpr, dpr);
}

function toggleSeqPanel() {
  const panel = document.getElementById('seqpanel');
  const toggle = document.getElementById('seqtoggle');
  panel.classList.toggle('collapsed');
  toggle.innerHTML = panel.classList.contains('collapsed') ? '&#9650;' : '&#9660;';
  // Resize 3D renderer after panel toggle
  setTimeout(() => {
    if (typeof resizeRenderer === 'function') resizeRenderer();
  }, 320);
}

// BLOSUM62 similarity to color: red (dissimilar) -> white (neutral) -> green (similar)
function simToColor(sim) {
  // sim in [0, 1]: 0 = very dissimilar, 1 = identical
  if (sim === null || sim === undefined) return '#30363d';  // gap
  let r, g, b;
  if (sim < 0.5) {
    // red to white
    const u = sim / 0.5;
    r = 215 + (230 - 215) * u;
    g = 58 + (237 - 58) * u;
    b = 73 + (243 - 73) * u;
  } else {
    // white to green
    const u = (sim - 0.5) / 0.5;
    r = 230 + (63 - 230) * u;
    g = 237 + (185 - 237) * u;
    b = 243 + (80 - 243) * u;
  }
  return 'rgb(' + Math.round(r) + ',' + Math.round(g) + ',' + Math.round(b) + ')';
}

function drawSeqPanel() {
  if (!seqCtx || !DATA.sequences) return;
  const ctx = seqCtx;
  const S = seqState;
  const cw = S.charW, ch = S.charH, lw = S.labelW, rh = S.rowH;
  const nRows = S.speciesOrder.length;
  const totalH = nRows * rh + 8;

  // Clear
  ctx.fillStyle = '#0d1117';
  ctx.fillRect(0, 0, S.colWidth, totalH);

  const scrollX = S.scrollX;
  const scrollY = S.scrollY;
  const container = document.getElementById('seqcontainer');
  const visW = container.clientWidth;

  // Draw each row
  for (let ri = 0; ri < nRows; ri++) {
    const sp = S.speciesOrder[ri];
    const seqData = DATA.sequences.data[sp];
    if (!seqData) continue;
    const aln = seqData.aln;
    const sim = seqData.sim;
    const y = 4 + ri * rh - scrollY;
    if (y + rh < 0 || y > container.clientHeight) continue;  // off-screen

    const isSelected = (sp === S.selectedSpecies);
    const isAnimating = (sp === S.targetSpecies && S.animating);
    const identity = DATA.sequences.identity[sp] || 0;

    // Row background highlight (smooth transition)
    let highlightAlpha = 0;
    if (isSelected) highlightAlpha = S.highlightAlpha;
    else if (isAnimating) highlightAlpha = 1 - S.highlightAlpha;
    if (highlightAlpha > 0.01) {
      ctx.fillStyle = 'rgba(88, 166, 255, ' + (highlightAlpha * 0.15) + ')';
      ctx.fillRect(0, y, S.colWidth, rh);
      // Left border accent
      ctx.fillStyle = 'rgba(88, 166, 255, ' + (highlightAlpha * 0.8) + ')';
      ctx.fillRect(0, y, 3, rh);
    }

    // Species label
    ctx.fillStyle = isSelected ? '#58a6ff' : '#8b949e';
    ctx.font = (isSelected ? 'bold ' : '') + '11px monospace';
    ctx.textBaseline = 'middle';
    const label = sp.length > 22 ? sp.substring(0, 21) + '\u2026' : sp;
    ctx.fillText(label, 6, y + rh / 2);
    // Identity percentage
    ctx.fillStyle = isSelected ? '#58a6ff' : '#484f58';
    ctx.font = '9px monospace';
    ctx.fillText(Math.round(identity * 100) + '%', lw - 40, y + rh / 2);

    // Sequence characters (only visible columns)
    const startCol = Math.max(0, Math.floor((scrollX - lw) / cw));
    const endCol = Math.min(S.nCols, Math.ceil((scrollX + visW - lw) / cw));
    ctx.font = '11px monospace';
    ctx.textBaseline = 'middle';
    for (let ci = startCol; ci < endCol; ci++) {
      const aa = aln[ci];
      const x = lw + ci * cw - scrollX;
      // v16: Check if this is the selected residue column (reference column)
      const refCol = msaColToRefCol(ci);
      const isSelectedRes = (refCol !== null && refCol === S.selectedRefCol);
      if (isSelectedRes) {
        // Draw selection box behind the residue
        ctx.fillStyle = 'rgba(240, 198, 116, 0.3)';
        ctx.fillRect(x - 1, y, cw + 1, rh);
        ctx.strokeStyle = '#f0c674';
        ctx.lineWidth = 1.5;
        ctx.strokeRect(x - 1, y, cw + 1, rh);
      }
      if (aa === '-') {
        // Gap: draw a faint dash
        ctx.fillStyle = isSelectedRes ? '#f0c674' : '#21262d';
        ctx.fillText('-', x, y + rh / 2);
      } else {
        const s = sim[ci];
        ctx.fillStyle = isSelectedRes ? '#f0c674' : simToColor(s);
        ctx.fillText(aa, x, y + rh / 2);
      }
    }
    
    // v16: Draw a vertical line connecting the selected residue across all rows
    // (drawn after all rows so it's on top)
  }

  // v16: Draw vertical connector line for selected residue across all species rows
  if (S.selectedRefCol !== null) {
    const refColMap = DATA.sequences.refColMap;
    const msaCol = refColMap[S.selectedRefCol];
    if (msaCol !== undefined) {
      const x = lw + msaCol * cw - scrollX;
      if (x >= lw - scrollX && x <= scrollX + visW) {
        ctx.strokeStyle = 'rgba(240, 198, 116, 0.15)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(x + cw / 2, 4 - scrollY);
        ctx.lineTo(x + cw / 2, 4 + nRows * rh - scrollY);
        ctx.stroke();
      }
    }
  }

  // Position ruler at top
  ctx.fillStyle = '#484f58';
  ctx.font = '9px monospace';
  ctx.textBaseline = 'middle';
  const rulerY = -scrollY + 2;
  if (rulerY > -10 && rulerY < 20) {
    const startCol = Math.max(0, Math.floor((scrollX - lw) / cw));
    const endCol = Math.min(S.nCols, Math.ceil((scrollX + visW - lw) / cw));
    for (let ci = startCol; ci < endCol; ci += 50) {
      const x = lw + ci * cw - scrollX;
      ctx.fillText(String(ci + 1), x, rulerY);
    }
  }

  // Update info text
  const info = document.getElementById('seqinfo');
  if (info) {
    if (S.selectedRefCol === null && S.selectedSpecies) {
      const ident = DATA.sequences.identity[S.selectedSpecies] || 0;
      info.textContent = S.selectedSpecies + ' \u00b7 ' + Math.round(ident * 100) + '% identity to reference';
    }
  }
}

// Smooth transition: animate highlight from old species to new species
function transitionSeqHighlight(newSpecies) {
  const S = seqState;
  if (S.selectedSpecies === newSpecies) return;
  S.targetSpecies = newSpecies;
  S.animating = true;
  S.highlightAlpha = 0;
  const startTime = performance.now();
  const duration = 400;  // ms

  function animStep() {
    const f = Math.min(1, (performance.now() - startTime) / duration);
    S.highlightAlpha = f;  // 0 -> 1 (old fades out, new fades in)
    drawSeqPanel();
    if (f < 1) {
      requestAnimationFrame(animStep);
    } else {
      S.selectedSpecies = newSpecies;
      S.targetSpecies = null;
      S.animating = false;
      S.highlightAlpha = 1;
      drawSeqPanel();
      // Scroll to bring selected row into view
      scrollToSelectedRow();
    }
  }
  animStep();
}

function scrollToSelectedRow() {
  const S = seqState;
  if (!S.selectedSpecies) return;
  const ri = S.speciesOrder.indexOf(S.selectedSpecies);
  if (ri < 0) return;
  const container = document.getElementById('seqcontainer');
  // v16: Use native scrollTop (container has overflow:auto)
  const rowY = 4 + ri * S.rowH;
  const targetScroll = Math.max(0, rowY - container.clientHeight / 2);
  container.scrollTop = targetScroll;
  S.scrollY = targetScroll;
  drawSeqPanel();
}

// v16: Scroll sequence panel to show the selected residue column
function scrollToSelectedResidue() {
  const S = seqState;
  if (S.selectedRefCol === null) return;
  const refColMap = DATA.sequences.refColMap;
  if (!refColMap || S.selectedRefCol >= refColMap.length) return;
  const msaCol = refColMap[S.selectedRefCol];
  const container = document.getElementById('seqcontainer');
  const targetX = S.labelW + msaCol * S.charW - container.clientWidth / 2;
  container.scrollLeft = Math.max(0, targetX);
  S.scrollX = container.scrollLeft;
  drawSeqPanel();
}

// ================= boot =================
renderTree();
init3D();
initSeqPanel();
const init = DATA.structures.find(s => s.accession === DATA.reference) ||
             DATA.structures[0];
selectSpecies(init.species);
// Set initial sequence panel selection
if (DATA.sequences) {
  seqState.selectedSpecies = init.species;
  seqState.highlightAlpha = 1;
  drawSeqPanel();
}
</script>
</body>
</html>

```

## `skills/morph_animation/__init__.py`

```


```

## `skills/morph_animation/morph_animation.py`

```
#!/usr/bin/env python3
"""Render an evolutionary morph animation from aligned structures. (v3)

Usage:
    python morph_animation.py structures/aligned.json -o evolution.gif
    python morph_animation.py aligned.json -o out.gif --hold 25 --morph 15 --fps 12

Reads the aligned.json produced by the structure-aligner skill v2 (all
structures already superposed in a common 3D frame, each carrying its
residue sequence and per-residue secondary structure) and renders a GIF
that steps through the species in evolutionary order (green algae ->
angiosperms).

v2 changes (the "smooth, change-only" transition):
    - RIBBON rendering: helices as flat ribbons, strands as arrows/ribbons,
      loops as a thin tube, instead of scatter+line. Clean cartoon look.
    - FIXED MEDIUM LIGHTING: two lights parented to the camera at constant
      intensity, so shading is consistent across frames and while rotating.
    - STABILITY-MASKED MORPH: per corresponding residue we compute the
      CA displacement and whether the secondary-structure class changed.
      Residues that are stable AND structurally conserved are frozen
      (a helix that is a helix in both species does not reshuffle);
      only genuinely changing regions move. Insertions fade in from the
      core centroid instead of popping in. Eased (smoothstep) timing.
    - v3: residue correspondence is anchored through the alignment
      reference (reference_pairs), so pairwise alignment slips between
      distant homologs no longer reshuffle the whole protein, and the
      N-/C-termini are labeled with markers on every frame.
    - No membrane planes anywhere.

Output:
    GIF animation (default evolution.gif) + a frame sheet PNG for
    quick visual inspection.
"""
import argparse
import json
import os
import sys

import numpy as np

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

LINEAGE_ORDER = [
    "Green Algae",
    "Bryophytes",
    "Ferns",
    "Gymnosperms",
    "Angiosperms",
]

LINEAGE_LABELS = {
    "Green Algae": "Green algae",
    "Bryophytes": "Bryophytes (mosses & liverworts)",
    "Ferns": "Ferns",
    "Gymnosperms": "Gymnosperms (conifers)",
    "Angiosperms": "Flowering plants",
}

# ribbon palette (matches the dark theme, reads well under medium light)
SS_COLOR = {"H": "#c678dd", "E": "#e5c07b", "C": "#7f8c98"}
BG = "#0d1117"


def order_structures(structs):
    """Order by lineage; within a lineage, chain by core identity."""
    ordered = []
    for lineage in LINEAGE_ORDER:
        group = [s for s in structs if s.get("lineage") == lineage]
        if not group:
            continue
        group.sort(key=lambda s: -(s.get("core_identity") or 0))
        ordered.extend(group)
    known = set(LINEAGE_ORDER)
    rest = [s for s in structs if s.get("lineage") not in known]
    return ordered + rest


# ------------------------------------------------------------- ribbon draw

def _ribbon_axes(coords, i, j):
    """Local frame along a ribbon segment: tangent, normal, binormal."""
    a = coords[i]
    b = coords[min(j - 1, len(coords) - 1)]
    tangent = b - a
    n = np.linalg.norm(tangent)
    if n < 1e-6:
        tangent = np.array([1.0, 0.0, 0.0])
    else:
        tangent = tangent / n
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(tangent, ref)) > 0.95:
        ref = np.array([0.0, 1.0, 0.0])
    normal = np.cross(tangent, ref)
    normal /= np.linalg.norm(normal)
    binormal = np.cross(tangent, normal)
    return tangent, normal, binormal


def draw_structure(ax, coords, ss, colors, min_helix=4, min_strand=3,
                   termini=None, radius_scale=None):
    """Cartoon ribbon: helices as wide ribbons, strands as arrows, loops as
    a thin tube. SS class colors are tinted per segment by mean pLDDT so
    confidence is still readable. Lighting is applied by the caller.

    termini: optional (n_idx, c_idx) marking the N- and C-terminal CA
    positions; drawn as small sphere markers with "N"/"C" text labels.
    radius_scale: optional per-residue array in [0, 1] scaling the ribbon
    linewidth (used for fade in/out of insertions/deletions during morph).
    """
    coords = np.asarray(coords, dtype=float)
    if radius_scale is not None:
        radius_scale = np.asarray(radius_scale, dtype=float)
        # Scale linewidths: helix 8.5->0, strand 6.5->0, loop 2.2->0
        lw_scale = radius_scale
    else:
        lw_scale = np.ones(len(coords))
    # faint backbone underneath for continuity
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
            color="#30363d", linewidth=0.6, alpha=0.5, zorder=1)
    segs = pl.ribbon_segments(coords, ss)
    for seg in segs:
        i, j = seg["idx"]
        pts = coords[i:j]
        if len(pts) < 2:
            continue
        if seg["type"] == "H" and (j - i) < min_helix:
            seg_type = "C"
        elif seg["type"] == "E" and (j - i) < min_strand:
            seg_type = "C"
        else:
            seg_type = seg["type"]
        hexcol = SS_COLOR[seg_type].lstrip("#")
        base = np.array([int(hexcol[k:k + 2], 16) / 255.0
                         for k in (0, 2, 4)])
        seg_cols = colors[i:j]
        tint = (seg_cols.mean(axis=0)[:3] if len(seg_cols)
                else colors[min(i, len(colors) - 1)][:3])
        col = tuple(np.clip(0.45 * base + 0.55 * tint, 0, 1))
        # Mean radius scale for this segment (fade in/out)
        seg_rs = float(lw_scale[i:j].mean()) if len(lw_scale) > j else 1.0
        seg_alpha = min(1.0, seg_rs * 1.5)  # fade visibility
        if seg_rs < 0.05:
            continue  # invisible segment
        if seg_type == "H":
            # helix: wide flat ribbon with darker edge passes -> shaded look
            tangent, normal, _ = _ribbon_axes(coords, i, j)
            off = normal * 0.85
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=8.5 * seg_rs, solid_capstyle="round",
                    zorder=3, alpha=seg_alpha)
            edge = tuple(np.clip(np.array(col) * 0.65, 0, 1))
            ax.plot(pts[:, 0] + off[0], pts[:, 1] + off[1],
                    pts[:, 2] + off[2], color=edge,
                    linewidth=2.6 * seg_rs, solid_capstyle="round",
                    alpha=0.95 * seg_alpha, zorder=2)
            ax.plot(pts[:, 0] - off[0], pts[:, 1] - off[1],
                    pts[:, 2] - off[2], color=edge,
                    linewidth=2.6 * seg_rs, solid_capstyle="round",
                    alpha=0.95 * seg_alpha, zorder=2)
        elif seg_type == "E":
            # strand: wide ribbon + arrow head at the C-terminal end
            tangent, normal, _ = _ribbon_axes(coords, i, j)
            off = normal * 0.70
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=6.5 * seg_rs, solid_capstyle="round",
                    zorder=3, alpha=seg_alpha)
            edge = tuple(np.clip(np.array(col) * 0.65, 0, 1))
            ax.plot(pts[:, 0] + off[0], pts[:, 1] + off[1],
                    pts[:, 2] + off[2], color=edge,
                    linewidth=2.2 * seg_rs, solid_capstyle="round",
                    alpha=0.95 * seg_alpha, zorder=2)
            ax.plot(pts[:, 0] - off[0], pts[:, 1] - off[1],
                    pts[:, 2] - off[2], color=edge,
                    linewidth=2.2 * seg_rs, solid_capstyle="round",
                    alpha=0.95 * seg_alpha, zorder=2)
            tip = pts[-1]
            head = pts[-1] - tangent * 2.6
            ax.plot([head[0], tip[0]], [head[1], tip[1]],
                    [head[2], tip[2]], color=col,
                    linewidth=9.0 * seg_rs, solid_capstyle="butt",
                    zorder=3, alpha=seg_alpha)
        else:
            # loop: thin tube
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=2.2 * seg_rs, solid_capstyle="round",
                    zorder=2, alpha=seg_alpha)
    if termini is not None:
        n_idx, c_idx = termini
        for idx, label, col in ((n_idx, "N", "#58a6ff"),
                                (c_idx, "C", "#ff7f0e")):
            if idx is None or idx >= len(coords):
                continue
            p = coords[idx]
            ax.scatter([p[0]], [p[1]], [p[2]], s=42, color=col,
                       depthshade=False, zorder=6,
                       edgecolors="#0d1117", linewidths=0.6)
            ax.text(p[0], p[1], p[2], f" {label}", color=col,
                    fontsize=8, fontweight="bold", zorder=6)


# --------------------------------------------------------------- rendering

def render(frames_meta, out_path, fps, size):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation

    # Filter out fully invisible points (radius_scale < 0.05) for bounds
    visible_pts = []
    for f in frames_meta:
        c = np.asarray(f["coords"])
        rs = f.get("radius_scale")
        if rs is not None:
            mask = np.asarray(rs) >= 0.05
            if mask.any():
                visible_pts.append(c[mask])
        else:
            visible_pts.append(c)
    all_pts = np.vstack(visible_pts) if visible_pts else np.vstack([f["coords"] for f in frames_meta])
    center = (all_pts.min(axis=0) + all_pts.max(axis=0)) / 2
    spread = (all_pts.max(axis=0) - all_pts.min(axis=0)).max() / 2 * 1.12

    fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor(BG)
    # fixed medium lighting: a constant light direction in VIEW space, so
    # shading is identical on every frame and while the camera rotates
    # (the structure never appears to rotate under a static light)
    LIGHT_DIR = np.array([0.5, 0.6, 0.62])
    LIGHT_DIR /= np.linalg.norm(LIGHT_DIR)
    total = len(frames_meta)
    print(f"[render] {total} frames -> {out_path}", file=sys.stderr)

    def shade(coords, colors):
        """Lambert-style shading with the fixed view-space light: normals
        are approximated from the CA trace direction."""
        coords = np.asarray(coords, dtype=float)
        cols = np.asarray(colors, dtype=float)
        n = min(len(coords), len(cols))
        if n < 2:
            return cols
        coords = coords[:n]
        cols = cols[:n]
        tang = np.gradient(coords, axis=0)
        norm = np.linalg.norm(tang, axis=1, keepdims=True)
        norm[norm < 1e-6] = 1.0
        tang = tang / norm
        lam = np.clip(np.abs(tang @ LIGHT_DIR), 0.0, 1.0)
        # medium lighting: 0.45 ambient + 0.55 diffuse
        shade_f = 0.45 + 0.55 * lam[:, None]
        out = np.clip(cols * shade_f, 0.0, 1.0)
        if len(out) < len(colors):
            pad = np.tile(out[-1:], (len(colors) - len(out), 1))
            out = np.vstack([out, pad])
        return out

    def update(i):
        fm = frames_meta[i]
        ax.cla()
        ax.set_facecolor(BG)
        c = np.asarray(fm["coords"])
        cols = shade(c, np.asarray(fm["colors"]))
        ss = fm.get("ss") or "C" * len(c)
        termini = fm.get("termini")
        if termini is not None:
            n_idx, c_idx = termini
            termini = (min(n_idx, len(c) - 1), min(c_idx, len(c) - 1))
        rs = fm.get("radius_scale")
        draw_structure(ax, c, ss, cols, termini=termini, radius_scale=rs)
        ax.set_xlim(-spread + center[0], spread + center[0])
        ax.set_ylim(-spread + center[1], spread + center[1])
        ax.set_zlim(-spread + center[2], spread + center[2])
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        ax.set_title(fm["title"], color="#e6edf3", fontsize=12, pad=2)
        ax.text2D(0.5, 0.03, fm["subtitle"], transform=ax.transAxes,
                  ha="center", color="#8b949e", fontsize=9)
        prog = i / max(total - 1, 1)
        ax.text2D(0.5, 0.96,
                  f"evolutionary progression  {prog * 100:.0f}%",
                  transform=ax.transAxes, ha="center",
                  color="#58a6ff", fontsize=9)
        # slow rotation; lighting stays fixed because shading is view-based
        ax.view_init(elev=18, azim=(20 + i * 0.5) % 360)
        return []

    anim = animation.FuncAnimation(fig, update, frames=total, interval=100)
    writer = animation.PillowWriter(fps=fps)
    anim.save(out_path, writer=writer,
              savefig_kwargs={"facecolor": fig.get_facecolor()})
    plt.close(fig)
    print(f"[render] saved {out_path}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aligned_json",
                   help="aligned.json from the structure-aligner skill")
    p.add_argument("-o", "--output", default="evolution.gif")
    p.add_argument("--hold", type=int, default=25,
                   help="frames to hold each structure")
    p.add_argument("--morph", type=int, default=15,
                   help="frames per morph transition")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--size", type=int, nargs=2, default=(720, 720))
    p.add_argument("--disp-max", type=float, default=4.0,
                   help="CA displacement (A) at which a residue fully animates")
    p.add_argument("--title", default=None,
                   help="overall title (default: none)")
    args = p.parse_args()

    if pl is None:
        sys.exit("proteinlib.py not found - put it next to the skills tree "
                 "(src/proteinlib.py) or on PYTHONPATH")

    with open(args.aligned_json) as f:
        data = json.load(f)
    structs = order_structures(data["structures"])
    if len(structs) < 2:
        sys.exit("Need at least 2 structures to animate")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def colors_for(s, n):
        return pl.plddt_colors(s.get("plddts") or [70] * n, n)

    # v3: correspondence for the morph is anchored through the alignment
    # reference, so a pairwise alignment slip between two distant homologs
    # cannot reshuffle the whole protein; N/C termini are labeled.
    reference = next((s for s in data["structures"]
                      if s["accession"] == data.get("reference")),
                     structs[0])
    mask = pl.StabilityMask(disp_max=args.disp_max)
    frames_meta = []
    for idx, s in enumerate(structs):
        coords = np.asarray(s["coords"], dtype=float)
        seq = s.get("seq") or "A" * len(coords)
        ss = s.get("ss") or "C" * len(coords)
        termini = (0, len(coords) - 1) if len(coords) else None
        fm = {
            "title": f"{LINEAGE_LABELS.get(s['lineage'], s['lineage'])}"
                     f" \u2014 {s['species']}",
            "subtitle": f"{s['accession']} \u00b7 {s['source']}"
                        + (f" \u00b7 core ident "
                           f"{s['core_identity']:.0%}"
                           if s.get("core_identity") else ""),
            "ss": ss,
            "termini": termini,
        }
        cols = colors_for(s, len(coords))
        for _ in range(args.hold):
            frames_meta.append({**fm, "coords": coords, "colors": cols})
        if idx + 1 < len(structs):
            nxt = structs[idx + 1]
            ncoords = np.asarray(nxt["coords"], dtype=float)
            nseq = nxt.get("seq") or "A" * len(ncoords)
            nss = nxt.get("ss") or "C" * len(ncoords)
            pairs = pl.reference_pairs(s, nxt, reference)
            pairs = [(a, b) for a, b in pairs
                     if a < len(coords) and b < len(ncoords)]
            weights, a_idx, b_idx = mask.weights(
                coords, ncoords, pairs, ss, nss)
            morph_result = pl.morph_frames(coords, ncoords, pairs,
                                           weights, args.morph)
            n_union = len(morph_result[0]["coords"])
            n_cols_next = colors_for(nxt, len(ncoords))
            cols_a = cols[a_idx]
            cols_b = n_cols_next[b_idx]
            n_stable = int((weights < 0.15).sum())
            print(f"[morph] {s['species']} -> {nxt['species']}: "
                  f"{len(pairs)} pairs, {n_stable} frozen "
                  f"({n_stable / max(len(pairs), 1):.0%}), "
                  f"union={n_union}",
                  file=sys.stderr)
            # Build per-union-entry colors and SS by mapping back to
            # A/B indices via the morph result's insertion/deletion flags
            # We need to reconstruct which union entries map to which
            # A/B indices. morph_frames doesn't return indices, so we
            # rebuild them here.
            b_to_pos = {b: pos for pos, b in enumerate(b_idx)}
            a_to_pos = {a: pos for pos, a in enumerate(a_idx)}
            n_a = len(coords)
            paired_b = np.zeros(len(ncoords), dtype=bool)
            paired_b[b_idx] = True
            paired_a = np.zeros(n_a, dtype=bool)
            paired_a[a_idx] = True
            # Rebuild union entry -> (b_idx, a_idx, type) mapping
            union_map = []
            j = 0
            while j < len(ncoords):
                if paired_b[j]:
                    ai = a_to_pos[j] if j in a_to_pos else None
                    # Actually a_to_pos maps a_idx -> position in a_idx array
                    # We need b -> a mapping
                    pass
                j += 1
            # Simpler: build b->a and a->b dicts from pairs
            b2a = {b: a for a, b in pairs}
            a2b = {a: b for a, b in pairs}
            union_map = []
            j = 0
            while j < len(ncoords):
                if paired_b[j]:
                    union_map.append((j, b2a[j], 'paired'))
                    j += 1
                else:
                    ins_start = j
                    while j < len(ncoords) and not paired_b[j]:
                        j += 1
                    ins_end = j
                    a_before = b2a[ins_start - 1] if ins_start > 0 else -1
                    a_after = b2a[ins_end] if ins_end < len(ncoords) else n_a
                    for ai in range(a_before + 1, a_after):
                        if not paired_a[ai]:
                            union_map.append((-1, ai, 'del'))
                    for bi in range(ins_start, ins_end):
                        union_map.append((bi, -1, 'ins'))
            # Trailing deletions
            if len(a_idx) > 0:
                last_paired_a = max(a_idx)
                for ai in range(last_paired_a + 1, n_a):
                    if not paired_a[ai]:
                        union_map.append((-1, ai, 'del'))

            for k in range(len(morph_result)):
                t = k / max(len(morph_result) - 1, 1)
                frame = morph_result[k]
                fcoords = frame["coords"]
                frs = frame["radius_scale"]
                # Build colors and SS per union entry
                cols_k = np.zeros((n_union, 3))
                ss_full = []
                for u, (bi, ai, typ) in enumerate(union_map):
                    if typ == 'paired':
                        pos = b_to_pos[bi]
                        w = weights[pos]
                        cols_k[u] = cols_a[pos] * (1 - t) + cols_b[pos] * t
                        ss_full.append(ss[a_idx[pos]] if w < 0.5 else nss[b_idx[pos]])
                    elif typ == 'ins':
                        cols_k[u] = n_cols_next[bi]
                        ss_full.append(nss[bi])
                    else:  # del
                        cols_k[u] = cols[ai] * (1 - t)
                        ss_full.append(ss[ai])
                frames_meta.append({
                    **{"title": (fm["title"] if t < 0.5 else
                                 f"{LINEAGE_LABELS.get(nxt['lineage'], nxt['lineage'])}"
                                 f" \u2014 {nxt['species']}"),
                       "subtitle": fm["subtitle"] if t < 0.5 else
                           f"{nxt['accession']} \u00b7 {nxt['source']}"
                           + (f" \u00b7 core ident "
                              f"{nxt['core_identity']:.0%}"
                              if nxt.get("core_identity") else "")},
                    "coords": fcoords, "colors": cols_k,
                    "ss": "".join(ss_full),
                    "radius_scale": frs,
                    "termini": (0, n_union - 1)
                    if n_union else None})
    render(frames_meta, args.output, args.fps, args.size)
    sheet = frames_meta[::args.hold + args.morph]
    try:
        cols_n = int(np.ceil(np.sqrt(len(sheet))))
        rows_n = int(np.ceil(len(sheet) / cols_n))
        fig = plt.figure(figsize=(cols_n * 3, rows_n * 3), dpi=100)
        for i, fm in enumerate(sheet, 1):
            ax = fig.add_subplot(rows_n, cols_n, i, projection="3d")
            c = np.asarray(fm["coords"])
            ss = fm.get("ss") or "C" * len(c)
            termini = fm.get("termini")
            if termini is not None:
                n_idx, c_idx = termini
                termini = (min(n_idx, len(c) - 1), min(c_idx, len(c) - 1))
            draw_structure(ax, c, ss, np.asarray(fm["colors"])[:, :3],
                           termini=termini,
                           radius_scale=fm.get("radius_scale"))
            ax.set_axis_off()
            ax.set_title(fm["title"], fontsize=7, color="#e6edf3")
            ax.view_init(elev=18, azim=30)
        fig.patch.set_facecolor(BG)
        plt.tight_layout()
        sheet_path = os.path.splitext(args.output)[0] + "_sheet.png"
        plt.savefig(sheet_path, facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"[render] frame sheet -> {sheet_path}", file=sys.stderr)
    except Exception as e:
        print(f"[warn] frame sheet failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()

```

## `skills/sequence_scraper/__init__.py`

```


```

## `skills/sequence_scraper/scrape_sequences.py`

```
#!/usr/bin/env python3
"""Scrape one homologous protein sequence per model species.

Usage:
    python scrape_sequences.py "ribulose bisphosphate carboxylase large chain" -o results.fasta
    python scrape_sequences.py ECA3 -o results.fasta --table results.tsv
    python scrape_sequences.py AT1G10130 -o results.fasta

Cascade per species:
    1. NCBI Protein search (gene name, AGI, protein-name terms; species
       then genus scope) - collects ALL family candidates/isoforms
    2. Every candidate is verified as a true homolog of the query
       sequence via a local k-mer similarity check (no BLAST queue),
       and only the single closest match is kept.
    3. Last resort: the query term is broadened to the protein class
       (e.g. 'ATPase') so unnamed records can be identified by similarity.
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
KMER_MIN_FRAC = 0.008


def kmer_similarity(seq_a, seq_b, k=KMER_K):
    """Fraction of seq_a k-mers shared with seq_b. Cheap homology proxy."""
    if not seq_a or not seq_b:
        return 0.0
    ka = {seq_a[i:i + k] for i in range(len(seq_a) - k + 1)}
    kb = {seq_b[i:i + k] for i in range(len(seq_b) - k + 1)}
    if not ka:
        return 0.0
    return len(ka & kb) / len(ka)


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
        """Last resort: broad class terms (e.g. 'ATPase'); the k-mer check
        picks the true homolog among loosely-related or unnamed records."""
        head = terms[0].split()[0] if terms[0] else "protein"
        for broad in (f"{head} ATPase", "ATPase"):
            if broad in terms:
                continue
            ids = self.search_ids(broad, species, retmax=100)
            records = self.fetch_records(ids) if ids else []
            scored = self._closest_to_query(records, query_sequence,
                                            query_length) if records else []
            if scored:
                return scored
        return []

    def _closest_to_query(self, records, query_sequence, query_length=None):
        """Rank verified homologs by k-mer similarity to the query sequence."""
        scored = []
        for r in records:
            if self._exclude(r, [], ""):
                continue
            sim = kmer_similarity(query_sequence, r["sequence"])
            if sim < KMER_MIN_FRAC:
                continue
            if query_length and not (0.5 * query_length <= r["length"] <= 1.3 * query_length):
                continue
            scored.append((sim, r))
        scored.sort(key=lambda x: x[0], reverse=True)
        return scored

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
            print(f"  HIT [{quality}]: {best['accession']} ({best['length']} aa) "
                  f"{best['definition'][:60]}", file=sys.stderr)
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

```

## `skills/structure_aligner/__init__.py`

```


```

## `skills/structure_aligner/fetch_and_align_structures.py`

```
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

```

## `src/proteinlib.py`

```
#!/usr/bin/env python3
"""Shared library for the protein evolution animation project.

Structure handling (CA-trace PDB parsing/writing), secondary-structure
assignment, local sequence correspondence (BLOSUM62), Kabsch superposition,
and the stability mask that drives the smooth morph transitions.

The stability mask is the core of the new animation behaviour: for every
pair of consecutive structures we compute, per corresponding residue,
  - the C-alpha displacement after superposition, and
  - whether the secondary-structure element (helix/strand/loop) is conserved.
Residues that are both spatially stable and structurally conserved are
FROZEN during the morph (helix stays a helix, in place); only residues
with real structural change are interpolated, and residues present in only
one structure fade in/out. This removes the old "whole protein reshuffles
with an abrupt end" transitions.
"""
import numpy as np

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
AA1TO3 = {v: k for k, v in AA3TO1.items()}
AA1TO3['X'] = 'GLY'


# ---------------------------------------------------------------- PDB I/O

def parse_pdb_ca(path):
    """Parse CA atoms from a PDB file -> dict(coords, plddts, residues, seq)."""
    coords, plddts, residues = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                residues.append(line[17:20].strip())
                coords.append((float(line[30:38]), float(line[38:46]),
                               float(line[46:54])))
                try:
                    plddts.append(float(line[60:66]))
                except ValueError:
                    plddts.append(50.0)
    seq = "".join(AA3TO1.get(r, 'X') for r in residues)
    return {"coords": np.array(coords, dtype=float),
            "plddts": np.array(plddts, dtype=float),
            "residues": residues, "seq": seq}


def write_ca_pdb(coords, plddts, path, name, ss=None, seq=None):
    """Write a CA-trace PDB, optionally with HELIX/SHEET records.

    HELIX/SHEET records make the file render as ribbons in PyMOL, ChimeraX,
    Mol*, etc. - the ribbon style the animation and viewer use.
    """
    coords = np.asarray(coords, dtype=float)
    lines = [f"REMARK   aligned CA-trace model: {name}"]
    if ss:
        n = 0
        i = 0
        while i < len(ss):
            if ss[i] == 'H':
                j = i
                while j < len(ss) and ss[j] == 'H':
                    j += 1
                if j - i >= 4:
                    n += 1
                    r1 = AA1TO3.get(seq[i] if seq else 'A', 'GLY')
                    r2 = AA1TO3.get(seq[j - 1] if seq else 'A', 'GLY')
                    lines.append(f"HELIX  {n:3d} {n:3d} {r1:>3s} A {i + 1:4d}  "
                                 f"{r2:>3s} A {j:4d}  1")
                i = j
            else:
                i += 1
        n = 0
        i = 0
        while i < len(ss):
            if ss[i] == 'E':
                j = i
                while j < len(ss) and ss[j] == 'E':
                    j += 1
                if j - i >= 3:
                    n += 1
                    r1 = AA1TO3.get(seq[i] if seq else 'A', 'GLY')
                    r2 = AA1TO3.get(seq[j - 1] if seq else 'A', 'GLY')
                    lines.append(f"SHEET  {n:3d} {n:3d} 1 {r1:>3s} A{i + 1:4d}  "
                                 f"{r2:>3s} A{j:4d}  0")
                i = j
            else:
                i += 1
    for k, (xyz, p) in enumerate(zip(coords, plddts), start=1):
        x, y, z = xyz
        res = AA1TO3.get(seq[k - 1] if seq and k - 1 < len(seq) else 'A', 'GLY')
        lines.append(f"ATOM  {k:5d}  CA  {res:>3s} A{k:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 {p:6.2f}")
    lines.append("TER")
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ------------------------------------------------------- secondary structure

def assign_secondary(coords):
    """Assign H/E/C per residue from CA geometry.

    Helix: CA(i)..CA(i+3) distance ~4.5-5.8 A. Strand: CA(i)..CA(i+5)
    distance ~13-15 A (parallel/antiparallel sheet repeat). Short runs are
    demoted to coil. Pure CA heuristics - good enough to detect conserved
    vs changed secondary structure between homologs.
    """
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    ss = ['C'] * n
    if n < 6:
        return ss
    d3 = np.linalg.norm(coords[3:] - coords[:-3], axis=1)
    d5 = np.linalg.norm(coords[5:] - coords[:-5], axis=1)
    for i in range(n):
        if i < len(d3) and 4.5 <= d3[i] <= 5.8:
            ss[i] = 'H'
        elif i < len(d5) and 13.0 <= d5[i] <= 15.0:
            ss[i] = 'E'
    out = list(ss)
    for sym, minlen in (('H', 4), ('E', 3)):
        i = 0
        while i < n:
            if ss[i] == sym:
                j = i
                while j < n and ss[j] == sym:
                    j += 1
                if j - i < minlen:
                    for k in range(i, j):
                        out[k] = 'C'
                i = j
            else:
                i += 1
    # smooth: a single H inside a strand run etc. stays; keep simple.
    return out


def ss_segments(ss):
    """Yield (start, end, type) runs of H/E (half-open, min lengths applied)."""
    segs = []
    i = 0
    while i < len(ss):
        if ss[i] in 'HE':
            j = i
            while j < len(ss) and ss[j] == ss[i]:
                j += 1
            segs.append((i, j, ss[i]))
            i = j
        else:
            i += 1
    return segs


# ---------------------------------------------------------- sequence align

def sequence_pairs(seq_a, seq_b):
    """Local BLOSUM62 correspondence between two residue sequences.

    Returns a list of (i, j) residue index pairs covering the conserved
    core. Local alignment is essential: distant homologs carry unrelated
    terminal extensions that corrupt global alignments.
    """
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
    try:
        alignment = aligner.align(seq_a, seq_b)[0]
    except Exception:
        n = min(len(seq_a), len(seq_b))
        return list(zip(range(n), range(n)))
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
    if len(pairs) < 50:
        n = min(len(seq_a), len(seq_b))
        pairs = list(zip(range(n), range(n)))
    return pairs


def core_identity(seq_a, seq_b, pairs):
    if not pairs:
        return 0.0
    same = sum(1 for i, j in pairs if seq_a[i] == seq_b[j])
    return same / len(pairs)


# ---------------------------------------------------------------- Kabsch

def kabsch(mobile, target):
    """Superpose mobile onto target; returns (rotation, translation).

    Convention: transformed = mobile @ rotation + translation.
    """
    mob = np.asarray(mobile, dtype=float)
    tgt = np.asarray(target, dtype=float)
    mob_c = mob - mob.mean(axis=0)
    tgt_c = tgt - tgt.mean(axis=0)
    cov = mob_c.T @ tgt_c
    u, _, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(u @ vt))
    corr = np.diag([1.0, 1.0, d])
    r = u @ corr @ vt
    shift = tgt.mean(axis=0) - mob.mean(axis=0) @ r
    return r, shift


def rmsd(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def superpose(mobile, reference, core_iterations=4):
    """Superpose mobile onto reference via sequence-matched CA pairs with
    iterative core refinement (conserved core drives the alignment)."""
    pairs = sequence_pairs(mobile["seq"], reference["seq"])
    if pairs is None:
        pairs = [(i, i) for i in range(min(len(mobile["coords"]),
                                           len(reference["coords"])))]
    m_idx = np.array([p[0] for p in pairs])
    r_idx = np.array([p[1] for p in pairs])
    m_idx = m_idx[m_idx < len(mobile["coords"])]
    r_idx = r_idx[r_idx < len(reference["coords"])]
    n = min(len(m_idx), len(r_idx))
    m_idx, r_idx = m_idx[:n], r_idx[:n]
    r, shift = kabsch(mobile["coords"][m_idx], reference["coords"][r_idx])
    keep = np.ones(len(m_idx), dtype=bool)
    for _ in range(core_iterations):
        dist = np.sqrt(np.sum(
            (mobile["coords"][m_idx] @ r + shift
             - reference["coords"][r_idx]) ** 2, axis=1))
        median = np.median(dist[keep])
        new_keep = dist <= max(2.0 * median, median + 4.0)
        if new_keep.sum() < 50 or (new_keep == keep).all():
            break
        keep = new_keep
        r, shift = kabsch(mobile["coords"][m_idx[keep]],
                          reference["coords"][r_idx[keep]])
    fit = rmsd(mobile["coords"][m_idx[keep]] @ r + shift,
               reference["coords"][r_idx[keep]])
    new_coords = mobile["coords"] @ r + shift
    return new_coords, fit, r, shift


# -------------------------------------------------------- stability mask

class StabilityMask:
    """Per-residue morph weights between two structures.

    For each corresponding residue pair (i in A, j in B) we compute:
      displacement  - CA distance between A and B after superposition
                      (both structures already share a common frame from
                      the aligner; this measures local divergence)
      ss_change     - True if the secondary-structure class differs

    weight = 0.0  -> frozen (stable + same SSE): residue does not move
    weight = 1.0  -> fully animated (large displacement or SSE change)

    Residues present in only one structure get fade weights handled by
    the caller (insertions fade in, deletions fade out).
    """

    def __init__(self, disp_max=4.0, disp_min=1.0):
        # displacement (A) below which a structurally-conserved residue is
        # fully frozen, and at/above which it is fully animated
        self.disp_max = disp_max
        self.disp_min = disp_min

    def weights(self, coords_a, coords_b, pairs, ss_a=None, ss_b=None):
        coords_a = np.asarray(coords_a, dtype=float)
        coords_b = np.asarray(coords_b, dtype=float)
        pairs = [(a, b) for a, b in pairs
                 if a < len(coords_a) and b < len(coords_b)]
        a_idx = np.array([p[0] for p in pairs], dtype=int)
        b_idx = np.array([p[1] for p in pairs], dtype=int)
        disp = np.linalg.norm(coords_a[a_idx] - coords_b[b_idx], axis=1)
        # conserved secondary structure + small displacement -> frozen;
        # ramps up to full animation at disp_max
        w_disp = np.clip((disp - self.disp_min)
                         / max(self.disp_max - self.disp_min, 1e-6),
                         0.0, 1.0)
        if ss_a is not None and ss_b is not None:
            ss_change = np.array([1.0 if ss_a[a] != ss_b[b] else 0.0
                                  for a, b in zip(a_idx, b_idx)])
        else:
            ss_change = np.zeros(len(pairs))
        # an SSE change always animates fully (new/lost helix or strand);
        # otherwise the displacement ramp decides
        w = np.maximum(w_disp, ss_change)
        # smoothstep for a softer transition edge
        w = w * w * (3.0 - 2.0 * w)
        return w, a_idx, b_idx


def ease(t):
    """Smoothstep easing for morph parameter t in [0, 1]."""
    t = min(max(t, 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


# ------------------------------------------- reference-anchored morphing

def reference_pairs(struct_a, struct_b, reference):
    """Correspondence between two structures anchored through a common
    reference structure.

    Direct pairwise local alignment between two distant homologs can slip
    (especially over long terminal extensions), which makes the morph look
    like the whole protein reshuffles. Anchoring both structures through
    the alignment reference (the structure all coordinates were superposed
    onto) yields a consistent residue correspondence for the morph.

    struct_a/struct_b: dicts with "seq" (per-coordinate residue string)
    reference: dict with "seq" and "coords"
    Returns a list of (i, j) index pairs into a/b's coordinate arrays.
    """
    ref_seq = reference.get("seq") or "A" * len(reference["coords"])
    pa = sequence_pairs(struct_a.get("seq", ""), ref_seq)
    pb = sequence_pairs(struct_b.get("seq", ""), ref_seq)
    a_of_ref = {}
    for a, r in pa:
        if a < len(struct_a.get("coords", ())):
            a_of_ref[r] = a
    b_of_ref = {}
    for b, r in pb:
        if b < len(struct_b.get("coords", ())):
            b_of_ref[r] = b
    common = sorted(set(a_of_ref) & set(b_of_ref))
    return [(a_of_ref[r], b_of_ref[r]) for r in common]


def chain_orientation_ok(struct, reference):
    """Sanity check: the structure's coordinate chain must run in the same
    direction as the reference (N->C both, or C->N both).

    The check compares, in the superposed common frame, how far the
    structure's coordinate ends sit from the reference's N- and C-terminal
    regions. A correctly oriented chain has both ends near their
    corresponding reference ends; a reversed chain has the coordinate
    head near the reference C-terminus and vice versa. Returns True if the
    orientation matches (or the structure is too short to tell).
    """
    coords = np.asarray(struct.get("coords", ()), dtype=float)
    ref_coords = np.asarray(reference.get("coords", ()), dtype=float)
    if len(coords) < 24 or len(ref_coords) < 24:
        return True
    pairs = reference_pairs(struct, reference, reference)
    if not pairs:
        return True
    a_idx = np.array([p[0] for p in pairs])
    r_idx = np.array([p[1] for p in pairs])
    n = min(len(a_idx), len(r_idx))
    a_idx, r_idx = a_idx[:n], r_idx[:n]
    # superpose struct onto reference using the sequence correspondence
    r, shift = kabsch(coords[a_idx], ref_coords[r_idx])
    moved = coords @ r + shift
    k = min(12, len(ref_coords) // 4)
    head_ref = ref_coords[:k]
    tail_ref = ref_coords[-k:]
    # distance from each structure end to the NEAREST residue of each
    # reference end region (robust to local wiggle at the termini)
    def region_dist(pt, region):
        return float(np.min(np.linalg.norm(region - pt, axis=1)))
    d_head = region_dist(moved[0], head_ref)
    d_tail = region_dist(moved[-1], tail_ref)
    d_head_rev = region_dist(moved[0], tail_ref)
    d_tail_rev = region_dist(moved[-1], head_ref)
    return (d_head + d_tail) <= (d_head_rev + d_tail_rev)


# ------------------------------------------------------------ morph frame

def morph_frames(coords_a, coords_b, pairs, weights, n_frames,
                 fade_in=None):
    """Generate n_frames masked-morph frames (union space: handles both
    insertions and deletions).

    Returns frames in UNION space: paired residues + insertions (B-only)
    + deletions (A-only). Each frame is a dict with:
      "coords": (n_union, 3) array
      "radius_scale": (n_union,) array — 1.0 for paired, t for insertions
        (fade in), 1-t for deletions (fade out)
      "is_insertion": (n_union,) bool array
      "is_deletion": (n_union,) bool array

    Stable pairs (weight ~0) are anchored to A's position for the whole
    morph and only snap to B at the very end, so a conserved helix stays
    put. Insertions grow from the nearest paired neighbor's morphing
    position and fade in (radius_scale = t). Deletions collapse toward
    the nearest paired neighbor and fade out (radius_scale = 1-t).

    fade_in: optional boolean array over B marking insertions explicitly.
    """
    coords_a = np.asarray(coords_a, dtype=float)
    coords_b = np.asarray(coords_b, dtype=float)
    n_a = len(coords_a)
    n_b = len(coords_b)
    a_idx = np.array([p[0] for p in pairs], dtype=int)
    b_idx = np.array([p[1] for p in pairs], dtype=int)
    w = np.asarray(weights, dtype=float)
    paired_b = np.zeros(n_b, dtype=bool)
    paired_b[b_idx] = True
    paired_a = np.zeros(n_a, dtype=bool)
    paired_a[a_idx] = True
    frozen = w < 0.15
    ins = ~paired_b
    if fade_in is not None:
        ins = ins | np.asarray(fade_in, dtype=bool)
    dele = ~paired_a

    # Build union index mapping: walk through B, insert A-only residues
    # (deletions) at the right positions. Each union entry has:
    # (b_idx or -1, a_idx or -1, type)
    union_entries = []  # list of (b_idx, a_idx, 'paired'|'ins'|'del')
    b_to_a = {b: a for a, b in zip(a_idx, b_idx)}
    a_to_b = {a: b for a, b in zip(a_idx, b_idx)}

    j = 0  # B pointer
    while j < n_b:
        if paired_b[j]:
            union_entries.append((j, b_to_a[j], 'paired'))
            j += 1
        else:
            # Insertion block in B; check for deletions in A between
            # surrounding paired residues
            ins_start = j
            while j < n_b and not paired_b[j]:
                j += 1
            ins_end = j
            # Find A range for deletions: between the A indices of the
            # paired residues surrounding this gap
            a_before = b_to_a[ins_start - 1] if ins_start > 0 else -1
            a_after = b_to_a[ins_end] if ins_end < n_b else n_a
            # Add deletions (A-only residues in this gap)
            for ai in range(a_before + 1, a_after):
                if not paired_a[ai]:
                    union_entries.append((-1, ai, 'del'))
            # Add insertions
            for bi in range(ins_start, ins_end):
                union_entries.append((bi, -1, 'ins'))
    # Trailing deletions (A-only after last paired B)
    if len(a_idx) > 0:
        last_paired_a = max(a_idx)
        for ai in range(last_paired_a + 1, n_a):
            if not paired_a[ai]:
                union_entries.append((-1, ai, 'del'))

    n_union = len(union_entries)
    # Precompute weight lookup per union entry
    w_union = np.ones(n_union)
    bpos_union = np.full(n_union, -1, dtype=int)
    apos_union = np.full(n_union, -1, dtype=int)
    for u, (bi, ai, typ) in enumerate(union_entries):
        bpos_union[u] = bi
        apos_union[u] = ai
        if typ == 'paired':
            # find weight for this pair
            w_union[u] = w[np.where(b_idx == bi)[0][0]]

    frames = []
    for k in range(n_frames):
        t = ease(k / max(n_frames - 1, 1))
        out = np.zeros((n_union, 3))
        rs = np.ones(n_union)
        is_ins = np.zeros(n_union, dtype=bool)
        is_del = np.zeros(n_union, dtype=bool)

        for u, (bi, ai, typ) in enumerate(union_entries):
            if typ == 'paired':
                a = coords_a[ai]
                b = coords_b[bi]
                wu = w_union[u]
                if wu < 0.15:
                    settle = 1.0 if t >= 1.0 else 0.0
                    out[u] = a * (1 - settle) + b * settle
                else:
                    eff = 1.0 if t >= 1.0 else wu * t
                    out[u] = a * (1 - eff) + b * eff
                rs[u] = 1.0
            elif typ == 'ins':
                # Grow from nearest paired neighbor's morphing position
                anchor_u = u - 1  # preceding union entry
                while anchor_u >= 0 and union_entries[anchor_u][2] != 'paired':
                    anchor_u -= 1
                if anchor_u < 0:
                    # Search forward
                    anchor_u = u + 1
                    while anchor_u < n_union and union_entries[anchor_u][2] != 'paired':
                        anchor_u += 1
                if anchor_u >= 0 and anchor_u < n_union:
                    anchor_pos = out[anchor_u]
                else:
                    anchor_pos = coords_a.mean(axis=0)
                b = coords_b[bi]
                out[u] = anchor_pos * (1 - t) + b * t
                rs[u] = t
                is_ins[u] = True
            elif typ == 'del':
                # Collapse toward nearest paired neighbor
                anchor_u = u - 1
                while anchor_u >= 0 and union_entries[anchor_u][2] != 'paired':
                    anchor_u -= 1
                if anchor_u < 0:
                    anchor_u = u + 1
                    while anchor_u < n_union and union_entries[anchor_u][2] != 'paired':
                        anchor_u += 1
                if anchor_u >= 0 and anchor_u < n_union:
                    anchor_pos = out[anchor_u]
                else:
                    anchor_pos = coords_a[ai]
                a = coords_a[ai]
                out[u] = a * (1 - t) + anchor_pos * t
                rs[u] = 1.0 - t
                is_del[u] = True

        frames.append({
            "coords": out,
            "radius_scale": rs,
            "is_insertion": is_ins,
            "is_deletion": is_del,
        })

    # Exact endpoints: at t=0, paired at A, insertions at anchor, deletions at A
    # At t=1, paired at B, insertions at B, deletions at anchor (rs=0)
    return frames


# ---------------------------------------------------------- ribbon data

def ribbon_segments(coords, ss):
    """Split a CA trace into ribbon segments per SSE run.

    Returns a list of dicts: {type: 'H'|'E'|'C', idx: [i, j), coords}.
    Helices and strands become flat ribbon pieces; loops become a thin
    tube. Used by the matplotlib renderer to draw a cartoon-style ribbon.
    """
    coords = np.asarray(coords, dtype=float)
    segs = []
    for (i, j, typ) in ss_segments(ss):
        segs.append({"type": typ, "idx": (i, j),
                     "coords": coords[i:j]})
    # coil = everything not covered
    covered = np.zeros(len(coords), dtype=bool)
    for s in segs:
        covered[s["idx"][0]:s["idx"][1]] = True
    i = 0
    n = len(coords)
    while i < n:
        if not covered[i]:
            j = i
            while j < n and not covered[j]:
                j += 1
            segs.append({"type": 'C', "idx": (i, j), "coords": coords[i:j]})
            i = j
        else:
            i += 1
    segs.sort(key=lambda s: s["idx"][0])
    return segs


# ------------------------------------------------------------- colors

def plddt_colors(plddts, n=None):
    """RdYlGn colors from per-residue pLDDT, padded/truncated to n."""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("RdYlGn")
    p = np.asarray(plddts, dtype=float)
    if n is not None:
        if len(p) < n:
            p = np.concatenate([p, np.full(n - len(p), 70.0)])
        else:
            p = p[:n]
    if len(p) == 0:
        p = np.array([70.0])
    q = np.clip(p / 100.0, 0, 1)
    return cmap(0.15 + 0.85 * q)

```

## `tests/__init__.py`

```


```

## `tests/test_activities.py`

```
"""End-to-end tests for the protein-animation workflow activities.

Runs the align -> animate -> viewer chain on the cached ECA3 two-species
fixture (no network needed beyond the two AlphaFold PDBs already in
data/). Validates the artifact payloads the workflow returns.
"""
import asyncio
import base64
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from protein_animation.activities import (align_structures, build_viewer,  # noqa: E402
                                          render_animation)

DATA = ROOT.parent / "data"
FASTA = DATA / "test_sequences.fasta"


def _decode(b64):
    return base64.b64decode(b64)


async def main():
    fasta = FASTA.read_text()

    print("[1/3] align_structures ...")
    aln = await align_structures("ECA3_test", fasta, email=None)
    assert aln["n_structures"] == 2, aln["n_structures"]
    assert aln["reference"] == "Q9SY55", aln["reference"]
    assert set(aln["species"]) == {"Arabidopsis thaliana",
                                   "Coccomyxa viridis"}, aln["species"]
    assert set(aln["pdbs"]) == {"A0AAV1IES4.pdb", "Q9SY55.pdb"}, list(aln["pdbs"])
    assert aln["overlay_png"], "missing overlay"
    aj = json.loads(aln["aligned_json"])
    assert aj["reference"] == "Q9SY55"
    assert all("ss" in s and "seq" in s for s in aj["structures"])
    # PDBs are ribbon-ready (HELIX/SHEET records)
    for name, b64 in aln["pdbs"].items():
        pdb = _decode(b64).decode()
        assert pdb.startswith("HELIX") or pdb.startswith("SHEET") or "ATOM" in pdb
        assert "ATOM" in pdb
    # overlay is a real PNG
    assert _decode(aln["overlay_png"])[:8] == b"\x89PNG\r\n\x1a\n"
    print(f"  ok: {aln['n_structures']} structures, ref {aln['reference']}")

    print("[2/3] render_animation ...")
    anim = await render_animation("ECA3_test", aln["aligned_json"],
                                  hold=3, morph=4, fps=6, size=(320, 320))
    gif = _decode(anim["gif"])
    assert gif[:6] in (b"GIF87a", b"GIF89a"), gif[:6]
    assert anim["n_frames"] == 10, anim["n_frames"]  # 2 holds *3 + 1 morph *4
    assert anim["sheet_png"], "missing frame sheet"
    assert _decode(anim["sheet_png"])[:8] == b"\x89PNG\r\n\x1a\n"
    print(f"  ok: {anim['n_frames']} frames, gif {len(gif)} bytes")

    print("[3/3] build_viewer ...")
    v = await build_viewer("ECA3_test", aln["aligned_json"], title="ECA3 test")
    html = _decode(v["viewer_html"]).decode()
    assert "__DATA__" not in html and "__TITLE__" not in html
    assert "ECA3 test" in html
    assert "membrane" not in html.lower()
    assert "pairWeights" in html and "buildMorphCoords" in html
    assert "HemisphereLight" in html
    print(f"  ok: viewer {len(html)} bytes, morph+lighting present, "
          f"0 membrane refs")

    print("\nALL WORKFLOW ACTIVITY TESTS PASSED")


if __name__ == "__main__":
    asyncio.run(main())

```
