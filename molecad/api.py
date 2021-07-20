from typing import List

from fastapi import FastAPI
from pydantic.json import pydantic_encoder

from .api_db import molprops_pipeline
from .api_types import Properties

app = FastAPI()


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Properties]:
    return molprops_pipeline(smiles, skip, limit)
