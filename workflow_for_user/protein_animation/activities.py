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
