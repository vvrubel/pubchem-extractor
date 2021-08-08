from typing import Any, Dict, List

from loguru import logger
from mongordkit.Search import substructure
from pydantic import NonNegativeInt, PositiveInt
from pymongo import MongoClient
from pymongo.cursor import Cursor
from rdkit import Chem

from molecad.errors import NoDatabaseRecordError
from molecad.settings import settings
from molecad.utils import timer

db = MongoClient(settings.mongo_uri)[settings.db_name]


def paging_pipeline(mol_lst: List[str], skip: int, limit: int) -> List[Dict[str, Any]]:
    match_ = {"$match": {"index": {"$in": mol_lst}}}
    project_ = {"$project": {"_id": 0, "index": 0}}
    skip_ = {"$skip": skip}
    limit_ = {"$limit": limit}
    return [match_, project_, skip_, limit_]


def summary_pipeline(mol_lst: List[str]) -> List[Dict[str, Any]]:
    match_ = {"$match": {"index": {"$in": mol_lst}}}
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
            "MolecularWeight": {
                "Average": {"$round": ["$AvgMolW", 2]},
                "StandardDeviation": {"$round": ["$StdMolW", 2]},
            },
            "XLogP": {
                "Average": {"$round": ["$AvgLogP", 2]},
                "StandardDeviation": {"$round": ["$StdLogP", 2]},
            },
            "HBondDonorCount": {
                "Average": {"$round": ["$AvgDonor", 2]},
                "StandardDeviation": {"$round": ["$StdDonor", 2]},
            },
            "HBondAcceptorCount": {
                "Average": {"$round": ["$AvgAcceptor", 2]},
                "StandardDeviation": {"$round": ["$StdAcceptor", 2]},
            },
            "RotatableBondCount": {
                "Average": {"$round": ["$AvgRotatable", 2]},
                "StandardDeviation": {"$round": ["$StdRotatable", 2]},
            },
            "AtomStereoCount": {
                "Average": {"$round": ["$AvgAtomStereo", 2]},
                "StandardDeviation": {"$round": ["$StdAtomStereo", 2]},
            },
            "BondStereoCount": {
                "Average": {"$round": ["$AvgBondStereo", 2]},
                "StandardDeviation": {"$round": ["$StdBondStereo", 2]},
            },
            "Volume3D": {
                "Average": {"$round": ["$AvgVol3D", 2]},
                "StandardDeviation": {"$round": ["$StdVol3D", 2]},
            },
        }
    }
    return [match_, group_, project_]


def run_search(smiles: str) -> List[str]:
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    search_results: List[str] = substructure.SubSearch(q_mol, db[settings.molecules])

    if len(search_results) == 0:
        raise NoDatabaseRecordError
    else:
        logger.success(f"Found {len(search_results)} compounds for requested pattern.")
        return search_results


@timer
def compound_search(smiles: str, skip: NonNegativeInt, limit: PositiveInt) -> List[Cursor]:
    result = run_search(smiles)
    pipeline = paging_pipeline(result, skip, limit)
    cursor = db[settings.properties].aggregate(pipeline)
    return list(cursor)


@timer
def compound_search_summary(smiles: str) -> List[Cursor]:
    result = run_search(smiles)
    pipeline = summary_pipeline(result)
    cursor = db[settings.properties].aggregate(pipeline)
    return list(cursor)
