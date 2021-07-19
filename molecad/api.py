from typing import List

from fastapi import FastAPI

from .api_db import simple_search, summary_search
from .api_types import Properties, Summary
from .settings import settings  # must be changed

app = FastAPI()

# settings = Settings(_env_file="prod.env", _env_file_encoding="utf-8")
db = settings.get_db()
properties, molecules, mfp_counts = settings.get_collections()
db_collections = properties, molecules, mfp_counts


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Properties]:
    return simple_search(smiles, skip, limit)


@app.get("/v1/compound/summary", response_model=List[Summary])
def get_compound_summary(smiles: str) -> List[Summary]:
    return summary_search(smiles)
