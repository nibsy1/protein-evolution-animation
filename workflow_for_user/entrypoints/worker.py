"""Worker entrypoint for the protein-animation workflow."""
import asyncio

import mistralai.workflows as workflows

from protein_animation.workflow import ProteinAnimationWorkflow


async def main():
    await workflows.run_worker([ProteinAnimationWorkflow])


if __name__ == "__main__":
    asyncio.run(main())
