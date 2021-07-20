from typing import Optional

from pydantic import BaseModel


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


class OutCompoundSummary(BaseModel):
    InputSmiles: str
    MolecularWeight: Statistics
    XLogP: Statistics
    HBondDonorCount: Statistics
    HBondAcceptorCount: Statistics
    RotatableBondCount: Statistics
    AtomStereoCount: Statistics
    BondStereoCount: Statistics
    Volume3D: Statistics
