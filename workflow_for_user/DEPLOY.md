# Deploying the protein-animation workflow (Workflows tab)

This folder is the complete `protein-animation` workflow project. You asked to
have it under the **Workflows** tab in your Mistral account. Here is exactly
how, step by step.

## What this is

A Mistral Workflow (Temporal-backed, durable) that takes one input —
`{"query": "TPC1"}` (gene name, AGI code like `AT4G03560`, or protein name) —
and returns the full artifact set: sequences, aligned PDBs, the evolutionary
morph animation GIF, a frame sheet, and an interactive HTML viewer.

```
scrape_sequences -> align_structures -> render_animation -> build_viewer
   (sequence-scraper)  (structure-aligner)  (morph-animation)  (interactive-viewer)
```

A supervisor durable agent writes a short markdown report at the end.

## Option A — Managed deployment (recommended): shows up under Workflows

Mistral managed deployments build this folder into a Docker image and run the
worker for you. The workflow then appears under **Workflows** in the Console
and you can click **Start Workflow** with `{"query": "TPC1"}`.

1. **Push this folder to a GitHub repo.**
   ```bash
   cd workflow_for_user
   git init
   git add .
   git commit -m "protein-animation workflow"
   # create a repo on github.com, then:
   git remote add origin https://github.com/<you>/protein-animation.git
   git push -u origin main
   ```
2. **Create the deployment in the Console.**
   - Open [console.mistral.ai](https://console.mistral.ai) → your workspace.
   - Sidebar: **Build → Workflows → Deployments → New deployment**.
   - Choose your GitHub repo + `main` branch.
   - **Build directory**: `.` (the repo root, since this folder *is* the project).
   - Dockerfile path: `Dockerfile` (default).
   - Create. Mistral Cloud builds the image and starts the worker. No API key
     in the repo — it authenticates with an injected service-account token.
3. **Run it.**
   - Sidebar: **Workflows**. You should see **protein-animation** ("Protein
     Animation").
   - Click it → **Start Workflow** → input `{"query": "TPC1"}`.
   - Watch it in the **Executions** tab. Total runtime ~9 min for a full
     16-species run.
4. **Optional**: bind a workspace secret if you want NCBI rate-limit headroom
   (not required). Add an `email` field to the input instead:
   `{"query": "TPC1", "email": "you@example.com"}`.

## Option B — Run a local worker (also registers under Workflows)

If you'd rather run the worker on your own machine, it registers the same
workflow in your workspace and it appears under **Workflows** in the Console
while the worker is running.

1. **Get a Mistral API key**: Console → your profile → **API keys** → create.
   (Shown once — copy it.)
2. **Install and run:**
   ```bash
   cd workflow_for_user
   pip install uv
   uv pip install --system -r requirements.txt   # or: python -m venv .venv && .venv/bin/pip install -r requirements.txt
   export MISTRAL_API_KEY=your-key-here          # NEVER commit this
   export DEPLOYMENT_NAME=protein-animation
   python entrypoints/worker.py                  # keep this terminal open
   ```
3. **Trigger from the Console**: Workflows → **protein-animation** →
   **Start Workflow** → `{"query": "TPC1"}`. The execution is routed to your
   local worker.
4. Stop the worker with `Ctrl+C` when done.

## Option C — Trigger via API (no Console)

```bash
curl -X POST https://api.mistral.ai/v1/workflows/protein-animation/execute \
  -H "Authorization: Bearer $MISTRAL_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"input": {"query": "TPC1"}}'
```

## Input reference

| Field | Type | Default | Description |
|---|---|---|---|
| `query` | string | — | **Required.** Gene name, AGI, or protein name |
| `email` | string | null | NCBI contact email (recommended) |
| `make_viewer` | bool | true | Also build the interactive HTML viewer |
| `hold` | int | 20 | Hold frames per species in the GIF |
| `morph` | int | 12 | Morph frames per transition |
| `fps` | int | 12 | GIF frames per second |

## Output

JSON with every artifact base64-encoded: `evolution_gif`,
`evolution_sheet_png`, `viewer_html`, `aligned_pdbs`, `overlay_png`,
`sequences_fasta`, `sequences_tsv`, `aligned_json`, `summary_tsv`, plus a
human-readable `report`.

## Project layout

```
protein_animation/
  activities.py     # 4 durable activities wrapping the skills
  workflow.py       # protein-animation workflow + supervisor agent
skills/             # the four skill scripts (imported as modules)
  sequence_scraper/ structure_aligner/ morph_animation/ interactive_viewer/
src/proteinlib.py   # shared geometry/morph library
entrypoints/worker.py
tests/test_activities.py
Dockerfile  requirements.txt  pyproject.toml  Makefile  README.md
```

## Tests

```bash
PYTHONPATH=. python tests/test_activities.py
```

Runs align → animate → viewer on the cached ECA3 fixture and validates every
artifact payload. (The fixture FASTA/PDBs are not in this folder; the test
expects them one level up under `data/`.)

## Why it's not already under Workflows

A workflow only appears under the **Workflows** tab once a **worker has
registered it** (Option B, needs your API key) or a **managed deployment is
serving it** (Option A, needs the GitHub repo). There is no "upload a folder"
button — the registration happens when the worker starts. This folder is
everything the worker needs; pick Option A or B above.
