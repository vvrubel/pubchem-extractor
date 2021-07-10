from typing import Any, Dict, List, Optional, Union

from fastapi import FastAPI
from pydantic import BaseModel

from molecad.data.core.db import Connect, db
from molecad.data.core.fingerprints import search

app = FastAPI()


class Molecule(BaseModel):
    CID: int
    MolecularFormula: str
    MolecularWeight: Union[float, str]  # почему-то в response приходит строка
    CanonicalSMILES: str
    InChI: str
    IUPACName: str
    XLogP: float
    HBondDonorCount: int
    HBondAcceptorCount: int
    RotatableBondCount: int
    AtomStereoCount: int
    BondStereoCount: int
    Volume3D: Optional[Union[float, int]]


@app.get("/v1/compound", response_model=List[Molecule])
def get_compounds(smiles: str, skip: int, limit: int) -> List[Molecule]:
    return search(smiles)[skip:][:limit]
