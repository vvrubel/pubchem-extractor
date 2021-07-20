from typing import Any, Dict, List, Optional

from mongordkit.Search import substructure
from pymongo.cursor import Cursor
from rdkit import Chem

from .errors import (
    BadRequestError,
    NoDatabaseRecordError,
)
from .settings import Settings


settings = Settings(_env_file="prod.env", _env_file_encoding="utf-8")
db = settings.get_db()
properties, molecules, mfp_counts = settings.get_collections()
db_collections = properties, molecules, mfp_counts


def paging_pipeline(skip: int, limit: int) -> List[Dict[str, int]]:
    skip_ = {"$skip": skip}
    limit_ = {"$limit": limit}
    return [skip_, limit_]


def summary_pipeline(smiles: str) -> List[Dict[str, Any]]:
    group_ = {
        "$group": {
            "_id": 0,
            "AvgMolW": {"$avg": "$MolecularWeight"},
            "StdMolW": {"$stdDevPop": "$MolecularWeight"},
            "AvgLogP": {"$avg": "$XLogP"},
            "StdLogP": {"$stdDevPop": "$XLogP"},
            "AvgDonor": {"$avg": "$HBondDonorCount"},
            "StdDonor": {"$stdDevPop": "$HBondDonorCount"},
            "AvgAcceptor": {"$avg": "$HBondAcceptorCount"},
            "StdAcceptor": {"$stdDevPop": "$HBondAcceptorCount"},
            "AvgRotatable": {"$avg": "$RotatableBondCount"},
            "StdRotatable": {"$stdDevPop": "$RotatableBondCount"},
            "AvgAtomStereo": {"$avg": "$AtomStereoCount"},
            "StdAtomStereo": {"$stdDevPop": "$AtomStereoCount"},
            "AvgBondStereo": {"$avg": "$BondStereoCount"},
            "StdBondStereo": {"$stdDevPop": "$BondStereoCount"},
            "AvgVol3D": {"$avg": "$Volume3D"},
            "StdVol3D": {"$stdDevPop": "$Volume3D"},
        }
    }
    project_ = {
        "$project": {
            "_id": 0,
            "InputSmiles": smiles,
            "MolecularWeight": {"Average": "$AvgMolW", "StandardDeviation": "$StdMolW"},
            "XLogP": {"Average": "$AvgLogP", "StandardDeviation": "$StdLogP"},
            "HBondDonorCount": {"Average": "$AvgDonor", "StandardDeviation": "$StdDonor"},
            "HBondAcceptorCount": {"Average": "$AvgAcceptor", "StandardDeviation": "$StdAcceptor"},
            "RotatableBondCount": {
                "Average": "$AvgRotatable",
                "StandardDeviation": "$StdRotatable",
            },
            "AtomStereoCount": {"Average": "$AvgAtomStereo", "StandardDeviation": "$StdAtomStereo"},
            "BondStereoCount": {"Average": "$AvgBondStereo", "StandardDeviation": "$StdBondStereo"},
            "Volume3D": {"Average": "$AvgVol3D", "StandardDeviation": "$StdVol3D"},
        }
    }
    return [group_, project_]


def search_route(*args) -> List[Dict[str, Any]]:
    route, smiles, skip, limit = args
    if route == "/v1/compound":
        return paging_pipeline(skip, limit)
    elif route == "/v1/compound/summary":
        return summary_pipeline(smiles)
    else:
        raise BadRequestError


def run_search(
    route: str, smiles: str, skip: Optional[int] = None, limit: Optional[int] = None
) -> List[Cursor]:
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    search_results: List[str] = substructure.SubSearch(q_mol, molecules)
    match_stage = [{"$match": {"index": {"$in": search_results}}}]
    routed_stages = search_route(route, smiles, skip, limit)
    pipeline = match_stage + routed_stages
    cursor: Cursor = properties.aggregate(pipeline)
    if not cursor:
        raise NoDatabaseRecordError
    else:
        return list(cursor)
