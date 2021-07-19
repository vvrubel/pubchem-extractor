from typing import Optional, Union

from pydantic import BaseModel


class Properties(BaseModel):
    CID: int
    MolecularFormula: str
    MolecularWeight: Union[str, float]  # TODO
    CanonicalSMILES: str
    InChI: str
    IUPACName: str
    XLogP: Optional[float] = None
    HBondDonorCount: Optional[int] = None
    HBondAcceptorCount: Optional[int] = None
    RotatableBondCount: Optional[int] = None
    AtomStereoCount: Optional[int] = None
    BondStereoCount: Optional[int] = None
    Volume3D: Optional[Union[float, int]]


# TODO я честно уберу все комментарии
class Summary(BaseModel):
    # Average_MolecularWeight: float
    # StandardDeviation_MolecularWeight: float
    Average_XLogP: float
    StandardDeviation_XLogP: float
    Average_HBondDonorCount: float
    StandardDeviation_HBondDonorCount: float
    Average_HBondAcceptorCount: float
    StandardDeviation_HBondAcceptorCount: float
    Average_RotatableBondCount: float
    StandardDeviation_RotatableBondCount: float
    Average_AtomStereoCount: float
    StandardDeviation_AtomStereoCount: float
    Average_BondStereoCount: float
    StandardDeviation_BondStereoCount: float
    Average_Volume3D: float
    StandardDeviation_Volume3D: float

# class Statistics(BaseModel):
#     Average: float
#     StandardDeviation: float


# class Summary(BaseModel):
#     MolecularWeight: List[Statistics]
#     XLogP: List[Statistics]
#     HBondDonorCount: List[Statistics]
#     HBondAcceptorCount: List[Statistics]
#     RotatableBondCount: List[Statistics]
#     AtomStereoCount: List[Statistics]
#     BondStereoCount: List[Statistics]
#     Volume3D: List[Statistics]
