



import json
import logging
from pathlib import Path

from maggma.stores import JSONStore, MemoryStore
from jobflow.core.store import JobStore

logging.basicConfig(level=logging.INFO)


store_path = Path("/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/data/endmembers/SER-Co/qhaflow/store.json")

store = JobStore(
    JSONStore(str(store_path), read_only=True),
    additional_stores={"data": MemoryStore()},
)

store.connect()

# Debug: list all jobs
logging.info("\nAll jobs in store:")
for doc in store.query(properties=["name", "uuid", "index"]):
    name = doc.get("name", "N/A")
    uuid = doc.get("uuid", "N/A")
    index = doc.get("index", "N/A")
    logging.info(f"  Name: {name} | Index: {index} | UUID: {uuid}")

inputs = ""
while inputs != "q":
    inputs = input("\nEnter UUID or q to quit: ").strip()
    if inputs == "q":
        break
    logging.info(f"Fetching output for UUID: {inputs}")
    try:
        job = store.query_one(
            criteria={"uuid": inputs},
            properties=["uuid", "index", "name"],
            sort={"index": -1},
        )
        if job is None:
            logging.warning(f"No job found with UUID: {inputs}")
            continue
        output = store.get_output(uuid=job["uuid"], which="last", load=True)
        # with open(f"output.json", "w", encoding="utf-8") as f:
        #     json.dump(output, f, indent=2)
        res_input = output["input"]
        with open("result.json", "w", encoding="utf-8") as f:
            json.dump(res_input, f, indent=2)
    except Exception as e:
        logging.error(f"Error fetching output for UUID {inputs}: {e}")
        continue

store.close()