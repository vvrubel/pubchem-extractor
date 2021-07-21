from typing import Any, Dict, List

from mongordkit.Search import substructure
from pymongo.cursor import Cursor
from rdkit import Chem

from .settings import settings
from .utils import timer

#settings = Settings(_env_file=".env.prod", _env_file_encoding="utf-8")
db = settings.get_db()
properties, molecules, mfp_counts = settings.get_collections()
db_collections = properties, molecules, mfp_counts


def paging_pipeline(skip: int, limit: int) -> List[Dict[str, int]]:
    skip_ = {"$skip": skip}
    limit_ = {"$limit": limit}
    return [skip_, limit_]


def summary_pipeline(smiles: str) -> List[Dict[str, Any]]:
    result = run_search(smiles)
    match_ = {"$match": {"index": {"$in": result}}}
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
    return [match_, group_, project_]


def run_search(smiles: str) -> List[str]:
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    search_results: List[str] = substructure.SubSearch(q_mol, molecules)
    return search_results


@timer
def compound_search(smiles_: str, skip_: int, limit_: int) -> List[Cursor]:
    result = run_search(smiles_)
    cursor = properties.find({"index": {"$in": result}}).skip(skip_).limit(limit_)
    return list(cursor)


@timer
def compound_search_summary(smiles: str) -> List[Cursor]:
    pipeline = summary_pipeline(smiles)
    cursor = properties.aggregate(pipeline)
    return list(cursor)
