from typing import List

from fastapi import FastAPI

from .api_db import
from .api_types import Properties, Summary

app = FastAPI()


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Properties]:
    return simple_search(smiles, skip, limit)


@app.get("/v1/compound/summary", response_model=List[Summary])
def get_compound_summary(smiles: str) -> List[Summary]:
    return summary_search(smiles)
