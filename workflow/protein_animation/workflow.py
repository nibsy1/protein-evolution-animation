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
