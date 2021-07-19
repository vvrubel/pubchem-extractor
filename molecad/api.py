from typing import List

from fastapi import FastAPI

from .api_db import simple_search
from .api_types import Properties

app = FastAPI()


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Properties]:
    return list(simple_search(smiles, skip, limit))
