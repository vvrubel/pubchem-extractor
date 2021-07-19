from typing import Optional

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
