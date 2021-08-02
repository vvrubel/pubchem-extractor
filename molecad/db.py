from typing import Any, Dict, List

from mongordkit.Search import similarity, substructure
from pydantic import NonNegativeInt, PositiveInt
from pymongo import MongoClient
from pymongo.cursor import Cursor
from rdkit import Chem
from rdkit.Chem import AllChem

from .settings import Settings
from .utils import timer

# Morgan fingerprints
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_LEN = 2048

# PyMongo configurations
DEFAULT_BATCH_SIZE = 100

setup = Settings()
db = MongoClient(setup.mongo_uri)[setup.db_name]
properties = db[setup.prop_collection]
molecules = db[setup.mol_collection]
mfp_counts = db[setup.mfp_collection]


def paging_pipeline(mol_lst: List[str], skip: int, limit: int) -> List[Dict[str, Any]]:
    """
    Функция принимает на вход список отфильтрованных молекул и параметры пагинации,
    подставляет эти значения в стадии пайплана и возвращает их список из функции.
    :param mol_lst: Список молекул из функции ``run_search``.
    :param skip: Число записей, которые нужно пропустить.
    :param limit: Число записей, которые нужно показать.
    :return: Список стадий.
    """
    match_ = {"$match": {"index": {"$in": mol_lst}}}
    project_ = {"$project": {"_id": 0, "index": 0}}
    skip_ = {"$skip": skip}
    limit_ = {"$limit": limit}
    return [match_, project_, skip_, limit_]


def summary_pipeline(mol_lst: List[str]) -> List[Dict[str, Any]]:
    """
    Функция принимает на вход список отфильтрованных молекул и строит пайплайн, в котором
    оставляет документы со smiles из входящего списка, затем группирует все
    документы и рассчитывает статистические параметры для числовых полей.
    :param mol_lst: Список молекул из функции ``run_search``.
    :return: Список стадий.
    """
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


def run_search(smiles: str):
    """
    Функция генерирует объект молекулы для работы rdkit, после этот объект используется для
    подструктурного поиска по коллекции "molecules".
    :param smiles: Строковое представление структуры молекулы.
    :return: Список молекул, удовлетворяют результатам поиска по заданной подструктуре.
    """
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    qfp = list(
        AllChem.GetMorganFingerprintAsBitVect(
            q_mol, DEFAULT_MORGAN_RADIUS, nBits=DEFAULT_MORGAN_LEN
        ).GetOnBits()
    )
    substructures: List[str] = substructure.SubSearch(q_mol, molecules)

    results = {}
    for doc in molecules.find({"index": {"$in": substructures}}):
        smi = doc["smiles"]
        tanimoto = similarity.calc_tanimoto(qfp, doc["fingerprints"]["morgan_fp"]["bits"])
        results[smi] = tanimoto

    return results


@timer
def compound_search(smiles: str, skip: NonNegativeInt, limit: PositiveInt) -> Cursor:
    q_mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    substructures = run_search(q_mol)

    pipeline = paging_pipeline(substructures, skip, limit)
    cursor = properties.aggregate(pipeline)
    return cursor


@timer
def compound_search_summary(smiles: str) -> Cursor:
    substructures = run_search(smiles)
    pipeline = summary_pipeline(substructures)
    cursor = properties.aggregate(pipeline)
    return cursor


if __name__ == "__main__":
    run_search("S(=O)(=O)NC(=O)N")
