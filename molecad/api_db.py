from typing import Any, Dict, Iterable, List, Tuple

from pymongo.cursor import Cursor
import rdkit.Chem
from pprint import pprint

from molecad.utils import chunked, url_encoder
from mongordkit.Search import similarity, substructure
from rdkit import Chem

from .settings import Settings, settings


# settings = Settings(_env_file="prod.env", _env_file_encoding="utf-8")
db = settings.get_db()
properties, molecules, mfp_counts = settings.get_collections()
db_collections = properties, molecules, mfp_counts


def run_search(smiles) -> List[str]:
    q_mol: rdkit.Chem.Mol = Chem.MolFromSmiles(smiles)
    return substructure.SubSearch(q_mol, molecules)


def simple_stages(smiles: str, skip: int, limit: int) -> List[Dict[str, Any]]:
    search_results = run_search(smiles)

    match_stage = {"$match": {"index": {"$in": search_results}}}
    skip_stage = {"$skip": skip}
    limit_stage = {"$limit": limit}
    return [match_stage, skip_stage, limit_stage]


def molprops_pipeline(smiles: str, skip: int, limit: int):
    pipeline = simple_stages(smiles, skip, limit)
    cursor = properties.aggregate(pipeline)
    return list(cursor)


if __name__ == "__main__":
    query = {"smiles": "CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N", "skip": 0, "limit": 2}
    molprops_pipeline("CC1=CC=C(C=C1)S(=O)(=O)NC(=O)N", 0, 2)
