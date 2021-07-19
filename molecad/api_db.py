from typing import Any, Dict, List

from mongordkit.Search import similarity, substructure
from rdkit import Chem

from .api import mfp_counts, molecules, properties
from .api_types import Properties, Summary


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


def simple_stages(result: List[str], skipped, limited) -> List[Dict[str, Any]]:
    match = {"$match": {"index": {"$in": result}}}
    skip = {"$skip": skipped}
    limit = {"$limit": limited}
    return [match, skip, limit]


def summary_stages(result: List[str]) -> List[Dict[str, Any]]:
    match = {"$match": {"index": {"$in": result}}}
    group = {
        "$group": {
            "_id": 0,
            "MolecularWeight.Average": {"$avg": "$MolecularWeight"},
            "MolecularWeight.StandardDeviation": {"$stdDevPop": "$MolecularWeight"},
            "XLogP.Average": {"$avg": "$XLogP"},
            "XLogP.StandardDeviation": {"$stdDevPop": "$XLogP"},
            "HBondDonorCount.Average": {"$avg": "$HBondDonorCount"},
            "HBondDonorCount.StandardDeviation": {"$stdDevPop": "$HBondDonorCount"},
            "HBondAcceptorCount.Average": {"$avg": "$HBondAcceptorCount"},
            "HBondAcceptorCount.StandardDeviation": {"$stdDevPop": "$HBondAcceptorCount"},
            "RotatableBondCount.Average": {"$avg": "$RotatableBondCount"},
            "RotatableBondCount.StandardDeviation": {"$stdDevPop": "$RotatableBondCount"},
            "AtomStereoCount.Average": {"$avg": "$AtomStereoCount"},
            "AtomStereoCount.StandardDeviation": {"$stdDevPop": "$AtomStereoCount"},
            "BondStereoCount.Average": {"$avg": "$BondStereoCount"},
            "BondStereoCount.StandardDeviation": {"$stdDevPop": "$BondStereoCount"},
            "Volume3D.Average": {"$avg": "$Volume3D"},
            "Volume3D.StandardDeviation": {"$stdDevPop": "$Volume3D"},
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


if __name__ == "__main__":
    smiles_ = "CCN1C=NC2=C(N=CN=C21)N"
    skip_ = 0
    limit_ = 10
