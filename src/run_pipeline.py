"""Full pipeline: taxonomy -> sequence search -> structures -> animation."""
from state.manager import WorkflowState
from agents.taxonomy_agent import TaxonomyAgent
from agents.search_agent import SearchAgent
from agents.structure_agent import StructureAgent
from agents.animation_agent import AnimationAgent
from config import Config


def main():
    print("=" * 60)
    print("\U0001f3ac PLANT PROTEIN EVOLUTION -> STRUCTURE ANIMATION")
    print("=" * 60)
    print(f"Target protein: {Config.TARGET_PROTEIN}")

    state = WorkflowState()

    if state.get_homologs_count() == 0:
        print("\n[1/3] Planning organisms and searching NCBI...")
        organisms = TaxonomyAgent().plan_organisms()
        if not organisms:
            print("\u274c Failed to plan organisms")
            return
        print(f"   Selected {len(organisms)} organisms")
        SearchAgent(state).execute_search(organisms)
        state.update_step(2)

    print(f"\n[2/3] Structures ({state.get_homologs_count()} homologs collected)...")
    structures = StructureAgent(state).run()
    if not structures:
        print("\u274c No structures obtained")
        return
    print(f"   {len(structures)} structures ready")

    print("\n[3/3] Rendering animation...")
    out = AnimationAgent(state).run()
    if out:
        print(f"\n\u2705 Animation saved: {out}")
    print("=" * 60)


if __name__ == "__main__":
    main()
