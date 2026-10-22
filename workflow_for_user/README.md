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
