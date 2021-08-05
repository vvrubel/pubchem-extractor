from typing import List, Optional

from pydantic import BaseModel, NonNegativeInt, PositiveInt


class Smiles(BaseModel):
    smiles: str


class PagingParams(BaseModel):
    skip: NonNegativeInt = 0
    limit: PositiveInt = 1


class SearchParams(BaseModel):
    smiles: str
    params: Optional[PagingParams]


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


class PageSearchModel(BaseModel):
    response: List[Compound]


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


class SummarySearchModel(BaseModel):
    response: List[CompoundSummary]
