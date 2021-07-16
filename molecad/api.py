from typing import List, Optional, Union

from fastapi import FastAPI
from pydantic import BaseModel

from molecad.db import sorting_search, stages, simple_search, summary_search

app = FastAPI()


class Properties(BaseModel):
    CID: int
    MolecularFormula: str
    MolecularWeight: Union[float, str]  # почему-то в response приходит строка
    CanonicalSMILES: str
    InChI: str
    IUPACName: str
    XLogP: Optional[float] = None
    HBondDonorCount: int
    HBondAcceptorCount: int
    RotatableBondCount: int
    AtomStereoCount: int
    BondStereoCount: int
    Volume3D: Optional[Union[float, int]]
    rdkit_index: str


@app.get("/v1/compound", response_model=List[Properties])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Properties]:
    return list(simple_search(smiles, skip, limit))



