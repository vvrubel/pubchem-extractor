import json
from typing import List

from fastapi import FastAPI
from pydantic.json import pydantic_encoder
from pymongo.cursor import Cursor

from .api_db import run_search
from .api_types import Properties, Summary

app = FastAPI()


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Cursor]:
    query = {"route": "/v1/compound", "smiles": smiles, "skip": skip, "limit": limit}
    res = run_search(**query)
    return res


@app.get("/v1/compound/summary", response_model=List[Summary])
def get_compound_summary(smiles: str) -> List[Cursor]:
    query = {
        "route": "/v1/compound/summary",
        "smiles": smiles,
    }
    res = run_search(**query)
    return res
