import json
import os
from datetime import datetime


class WorkflowState:
    """Persistent workflow state stored in state.json."""

    def __init__(self, path="state.json"):
        self.path = path
        if os.path.exists(self.path):
            with open(self.path) as f:
                self.data = json.load(f)
        else:
            self.data = {
                "step": 1,
                "created": datetime.now().isoformat(),
                "organisms": [],
                "homologs": [],
                "structures": [],
            }
            self.save()

    def save(self):
        with open(self.path, "w") as f:
            json.dump(self.data, f, indent=2)

    def update_step(self, step):
        self.data["step"] = step
        self.save()

    def add_organisms(self, organisms):
        self.data["organisms"] = organisms
        self.save()

    def add_homologs(self, homologs):
        self.data["homologs"].extend(homologs)
        self.save()

    def add_structures(self, structures):
        self.data["structures"] = structures
        self.save()

    def get_homologs_count(self):
        return len(self.data.get("homologs", []))

    def get_homologs(self):
        return self.data.get("homologs", [])

    def get_organisms(self):
        return self.data.get("organisms", [])
