from typing import List

from fastapi import FastAPI, Query
from pymongo.cursor import Cursor

from .api_db import compound_search, compound_search_summary
from .api_types import Compound, OutCompoundSummary

app = FastAPI()


@app.get("/v1/compound", response_model=List[Compound])
def get_compounds(
    smiles: str = Query(...), skip: int = Query(..., ge=0), limit: int = Query(..., ge=1)
) -> List[Cursor]:
    res = compound_search(smiles, skip, limit)
    return res


@app.get("/v1/compound/summary", response_model=List[OutCompoundSummary])
def get_compound_summary(smiles: str = Query(...)) -> List[Cursor]:
    res = compound_search_summary(smiles)
    return res
