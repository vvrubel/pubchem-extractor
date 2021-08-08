from typing import Any, Dict, List

from loguru import logger
from mongordkit.Search import substructure
from pydantic import NonNegativeInt, PositiveInt

# TODO: change to motor
from pymongo.cursor import Cursor
from rdkit import Chem

from .errors import NoDatabaseRecordError
from .settings import settings
from .utils import timer

# TODO: startup on event
db = settings.mongo_client_sync()


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
            "Found": {"$sum": 1},
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
            "Found": 1,
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


def search_substructures(smiles: str) -> List[str]:
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    search_results: List[str] = substructure.SubSearch(q_mol, db[settings.molecules])

    if len(search_results) == 0:
        raise NoDatabaseRecordError
    else:
        logger.success(f"Found {len(search_results)} compounds for requested pattern.")
        return search_results


# TODO: sorted results
@timer
def run_page_search(smiles: str, skip: NonNegativeInt, limit: PositiveInt) -> List[Cursor]:
    logger.info(f"Searching molecules for the substructure: {smiles}")
    mol_lst = search_substructures(smiles)
    logger.info(f"Applying the pagination parameters: skip – {skip}, limit – {limit}")
    pipeline = paging_pipeline(mol_lst, skip, limit)
    cursor = db[settings.properties].aggregate(pipeline)
    return list(cursor)


@timer
def run_summary_search(smiles: str) -> List[Cursor]:
    logger.info(f"Searching molecules for substructure: {smiles}")
    mol_lst = search_substructures(smiles)
    logger.info("Calculating the statistics for searched molecules.")
    pipeline = summary_pipeline(mol_lst)
    cursor = db[settings.properties].aggregate(pipeline)
    return list(cursor)
