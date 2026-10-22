import os
from dotenv import load_dotenv

load_dotenv()

class Config:
    MISTRAL_MODEL = "mistral-large-latest"
    NCBI_EMAIL = os.getenv("NCBI_EMAIL")
    TARGET_PROTEIN = os.getenv("TARGET_PROTEIN", "RuBisCO")
    
    # Evolutionary lineages to cover
    EVOLUTIONARY_STAGES = [
        "Chlorophyta",      # Green Algae (simple)
        "Bryophyta",        # Mosses
        "Pteridophyta",     # Ferns
        "Pinophyta",        # Gymnosperms
        "Magnoliophyta"     # Flowering Plants (complex)
    ]
    
    # How many species per lineage
    SPECIES_PER_LINEAGE = 2
