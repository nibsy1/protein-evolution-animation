PYTHON ?= python

.PHONY: install start-worker execute

install:
	uv pip install --system -r requirements.txt

start-worker:
	PYTHONPATH=. $(PYTHON) entrypoints/worker.py

execute:
	$(PYTHON) -c "import json,sys; from mistralai import Mistral; import os; \
	c = Mistral(api_key=os.environ['MISTRAL_API_KEY']); \
	e = c.workflows.execute_workflow(workflow_identifier='$(workflow)', input=json.loads('$(input)')); \
	print(e.model_dump_json(indent=2))"
