from typing import List

from fastapi import APIRouter
from pydantic import NonNegativeInt, PositiveInt

from .main import app
from .search import compound_search, compound_search_summary

router = APIRouter()


@app.get("/compound", response_model=List[Compound])
def get_compounds(smiles: str, skip: NonNegativeInt, limit: PositiveInt):
    res = compound_search(smiles, skip, limit)
    return list(res)


@app.get("/compound/summary", response_model=List[CompoundSummary])
def get_compound_summary(smiles: str):
    res = compound_search_summary(smiles)
    return list(res)
