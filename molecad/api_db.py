from typing import Any, Dict, List

from mongordkit.Search import similarity, substructure
from rdkit import Chem

from .api_types import Properties, Summary
from .settings import settings  # TODO must be changed to Settings

# TODO
# settings = Settings(_env_file="prod.env", _env_file_encoding="utf-8")
db = settings.get_db()
properties, molecules, mfp_counts = settings.get_collections()
db_collections = properties, molecules, mfp_counts


def sorted_results(smiles: str) -> List[str]:
    q_mol = Chem.MolFromSmiles(smiles)
    substr_set = set(substructure.SubSearch(q_mol, molecules, chirality=False))
    similar_search = similarity.SimSearchAggregate(q_mol, molecules, mfp_counts, 0.8)
    res = []
    for score, smi in sorted(similar_search, reverse=True):
        res.append(smi)
        if smi in substr_set:
            substr_set.remove(smi)
    res.extend(substr_set)
    return res


def simple_stages(result: List[str], skipped: int, limited: int) -> List[Dict[str, Any]]:
    match = {"$match": {"index": {"$in": result}}}
    skip = {"$skip": skipped}
    limit = {"$limit": limited}
    return [match, skip, limit]


def summary_stages(result: List[str]) -> List[Dict[str, Any]]:
    match = {"$match": {"index": {"$in": result}}}
    group = {
        "$group": {
            "_id": 0,
            # TODO
            # "Average_MolecularWeight": {"$avg": "$MolecularWeight"},
            # "StandardDeviation_MolecularWeight": {"$stdDevPop": "$MolecularWeight"},
            "Average_XLogP": {"$avg": "$XLogP"},
            "StandardDeviation_XLogP": {"$stdDevPop": "$XLogP"},
            "Average_HBondDonorCount": {"$avg": "$HBondDonorCount"},
            "StandardDeviation_HBondDonorCount": {"$stdDevPop": "$HBondDonorCount"},
            "Average_HBondAcceptorCount": {"$avg": "$HBondAcceptorCount"},
            "StandardDeviation_HBondAcceptorCount": {"$stdDevPop": "$HBondAcceptorCount"},
            "Average_RotatableBondCount": {"$avg": "$RotatableBondCount"},
            "StandardDeviation_RotatableBondCount": {"$stdDevPop": "$RotatableBondCount"},
            "Average_AtomStereoCount": {"$avg": "$AtomStereoCount"},
            "StandardDeviation_AtomStereoCount": {"$stdDevPop": "$AtomStereoCount"},
            "Average_BondStereoCount": {"$avg": "$BondStereoCount"},
            "StandardDeviation_BondStereoCount": {"$stdDevPop": "$BondStereoCount"},
            "Average_Volume3D": {"$avg": "$Volume3D"},
            "StandardDeviation_Volume3D": {"$stdDevPop": "$Volume3D"},
        }
    }
    return [match, group]


def simple_search(smiles: str, skip: int, limit: int) -> List[Properties]:
    res = sorted_results(smiles)
    pipeline = simple_stages(res, skip, limit)
    cursor = properties.aggregate(pipeline)
    return list(cursor)


def summary_search(smiles) -> List[Summary]:
    res = sorted_results(smiles)
    pipeline = summary_stages(res)
    cursor = properties.aggregate(pipeline)
    return list(cursor)
