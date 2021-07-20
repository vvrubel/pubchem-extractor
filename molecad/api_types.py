from typing import List, Optional

from pydantic import BaseModel


class Properties(BaseModel):
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
    Average: float
    StandardDeviation: float


class Summary(BaseModel):
    MolecularWeight: List[Statistics]
    XLogP: List[Statistics] = []
    HBondDonorCount: List[Statistics]
    HBondAcceptorCount: List[Statistics]
    RotatableBondCount: List[Statistics]
    AtomStereoCount: List[Statistics]
    BondStereoCount: List[Statistics]
    Volume3D: List[Statistics] = []
