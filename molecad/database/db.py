from typing import Any, Dict, List, Tuple

from loguru import logger
from mongordkit.Database import registration
from pymongo.collection import Collection
from pymongo.database import Database
from rdkit import Chem

from molecad.cli.errors import EmptySmilesError


def drop_db(db: Database) -> None:
    lst = db.list_collection_names()
    if len(lst) > 0:
        for item in lst:
            db.drop_collection(item)
            logger.info(f"Коллекция {item} удалена.")


def create_indexes(*args: Collection) -> None:
    for arg in args:
        arg.create_index("CID", unique=True)
        logger.info(f'На коллекции {arg.name} создан уникальный индекс "CID".')
        arg.create_index("index")
        logger.info(f'На коллекции {arg.name} создан индекс – "index".')


def register_from_smiles(smiles: str) -> Dict[str, Any]:
    if smiles is None:
        raise EmptySmilesError

    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise TypeError

    scheme = registration.MolDocScheme()
    return scheme.generate_mol_doc(rdmol)


def create_molecules(
    data: List[Dict[str, Any]], mol_collection: Collection
) -> Tuple[List[Dict[str, Any]], int]:
    inserted = 0
    for molprop in data:
        if "CanonicalSMILES" not in molprop:
            continue

        smiles = molprop["CanonicalSMILES"]
        cid = molprop["CID"]

        try:
            scheme = register_from_smiles(smiles)
        except (TypeError, EmptySmilesError):
            continue
        else:
            # добавляем в схему поле "CID"
            scheme["CID"] = cid

        mol_collection.insert_one(scheme)
        inserted += 1
        # к исходным данным добавляем поле "index" из схемы
        molprop["index"] = scheme["index"]

        # приводим значение ключа "MolecularWeight" к float
        mol_weight = molprop["MolecularWeight"]
        if mol_weight is not None:
            molprop["MolecularWeight"] = float(mol_weight)

    return data, inserted


def upload_data(data: List[Dict[str, Any]], prop_collection: Collection) -> Tuple[int, int]:
    res = prop_collection.insert_many(data, ordered=False).inserted_ids
    return len(res), (len(data) - len(res))


def delete_broken(prop_collection: Collection) -> Tuple[int, int]:
    without_smiles = prop_collection.delete_many(
        {"CanonicalSMILES": {"$exists": False}}
    ).deleted_count
    without_schemes = prop_collection.delete_many({"index": {"$exists": False}}).deleted_count
    return without_smiles, without_schemes
