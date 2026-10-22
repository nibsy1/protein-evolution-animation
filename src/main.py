from state.manager import WorkflowState
from agents.taxonomy_agent import TaxonomyAgent
from agents.search_agent import SearchAgent
from config import Config

def main():
    print("=" * 60)
    print("\U0001f331 PLANT PROTEIN EVOLUTION WORKFLOW")
    print("=" * 60)
    print(f"Target Protein: {Config.TARGET_PROTEIN}")
    print(f"Model: {Config.MISTRAL_MODEL}")
    print("=" * 60)
    
    # Initialize state
    state = WorkflowState()
    print(f"\n\U0001f4c1 Current workflow step: {state.data['step']}")
    
    # Check if Step 1 already done
    if state.data['step'] > 1:
        print("\u26a0 Step 1 already completed!")
        print(f"   Found {state.get_homologs_count()} homologs in state.json")
        choice = input("   Run anyway? (y/n): ")
        if choice.lower() != 'y':
            return
    
    # STEP 1: Plan organisms
    print("\n" + "=" * 60)
    print("STEP 1: Planning Evolutionary Coverage")
    print("=" * 60)
    
    tax_agent = TaxonomyAgent()
    organisms = tax_agent.plan_organisms()
    
    if not organisms:
        print("\u274c Failed to get organism list. Check API key.")
        return
    
    print(f"\n\U0001f4cb Selected {len(organisms)} organisms:")
    for org in organisms:
        print(f"   \u2022 {org['name']} ({org['lineage']})")
    
    # STEP 2: Search for sequences
    print("\n" + "=" * 60)
    print("STEP 2: Searching Protein Databases")
    print("=" * 60)
    
    search_agent = SearchAgent(state)
    search_agent.execute_search(organisms)
    
    # Mark step complete
    state.update_step(2)
    
    # Summary
    print("\n" + "=" * 60)
    print("\u2705 STEP 1 COMPLETE!")
    print("=" * 60)
    print(f"Total homologs collected: {state.get_homologs_count()}")
    print(f"Next step: AlphaFold prediction (Step 2)")
    print(f"Data saved to: {state.path}")
    print("=" * 60)

if __name__ == "__main__":
    main()
