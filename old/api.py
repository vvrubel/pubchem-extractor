from typing import List, Optional

from fastapi import FastAPI
from pydantic import BaseModel, NonNegativeInt, PositiveInt

from .api_db import compound_search, compound_search_summary
from src.settings import settings


class Compound(BaseModel):
    CID: int
    MolecularFormula: str
    MolecularWeight: float
    CanonicalSMILES: str
    InChI: str
    IUPACName: Optional[str] = None
    XLogP: Optional[float] = None
    HBondDonorCount: int
    HBondAcceptorCount: int
    RotatableBondCount: int
    AtomStereoCount: int
    BondStereoCount: int
    Volume3D: Optional[float] = None


class Statistics(BaseModel):
    Average: Optional[float] = None
    StandardDeviation: Optional[float] = None


class CompoundSummary(BaseModel):
    MolecularWeight: Statistics
    XLogP: Statistics
    HBondDonorCount: Statistics
    HBondAcceptorCount: Statistics
    RotatableBondCount: Statistics
    AtomStereoCount: Statistics
    BondStereoCount: Statistics
    Volume3D: Statistics


app = FastAPI()


@app.get(f"/{settings.api_version}/compound", response_model=List[Compound])
def get_compounds(smiles: str, skip: NonNegativeInt, limit: PositiveInt):
    res = compound_search(smiles, skip, limit)
    return list(res)


@app.get(f"/{settings.api_version}/compound/summary", response_model=List[CompoundSummary])
def get_compound_summary(smiles: str):
    res = compound_search_summary(smiles)
    return list(res)
