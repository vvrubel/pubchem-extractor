from typing import Any, Dict, List, Tuple

import pymongo
from loguru import logger
from mongordkit.Database import registration
from mongordkit.Search import similarity, substructure
from pymongo.collection import Collection
from pymongo.database import Database
from rdkit import Chem

from molecad.settings import settings


def prepare_db(db: Database, drop: bool) -> None:
    """
    Если в базе данных есть коллекции, то она будет очищена, после чего коллекции будут созданы
    заново и на них будет создан индекс ``CID``.
    :param db: Объект базы данных.
    :param drop: Флаг для очистки базы данных.
    :return: None.
    """
    lst = db.list_collection_names()
    if drop and len(lst) > 0:
        for item in lst:
            db.drop_collection(item)
            logger.info(f"Коллекция {item} удалена.", fg="yellow")

    names = [settings.properties, settings.molecules]
    for name in names:
        collection: pymongo.collection.Collection = db.create_collection(name)
        collection.create_index("CID", unique=True)
        logger.info(f'На коллекции {collection.name} создан индекс "CID".')


def register_from_smiles(smiles: str) -> Dict[str, Any]:
    """
    Создает объект схемы, определяющий модель документа молекул, и назначает индекс на коллекции.
    :param smiles: Строковое представление молекулы.
    :return: Схема документа молекулы.
    """
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise TypeError
    scheme = registration.MolDocScheme()
    scheme.set_index("CanonicalSmiles")
    return scheme.generate_mol_doc(rdmol)


def populate_w_schema(
    data: List[Dict[str, Any]], mol_collection: Collection
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Функция вызывается внутри команды ``populate`` и получает список из словарей (в будущем
    документов): в нем смотрит на наличие ключа ``CanonicalSMILES`` и, если он имеется, генерирует
    схему документа молекулы, добавляя в схему ключ ``CID``, вставляет документ в базу в случае
    успешного завершения предыдущих этапов. В последнюю очередь создает на коллекции ``properties``
    индекс ``rdkit_index``.
    :param data: Данные, загруженные из файла ``.json``.
    :param mol_collection: Коллекция, в которую будут вставлены документы, являющиеся
    представлением молекул.
    :return: Кортеж из прешедших в функцию данных, к которым добавлено поле ``CID``, и количество
    сгенерированных схем.
    """
    n = 0
    for mol in data:
        if "CanonicalSMILES" not in mol:
            continue
        smiles = mol["CanonicalSMILES"]
        cid = mol["CID"]
        try:
            scheme = register_from_smiles(smiles)
        except TypeError:
            # if scheme is None:
            continue
        scheme["CID"] = cid

        mol_collection.insert_one(scheme)
        n += 1

        mol["rdkit_index"] = scheme["index"]
    return data, n


def upload_data(
    data: List[Dict[str, Any]], collection: pymongo.collection.Collection
) -> Tuple[int, int]:
    """
    Загружает данные в коллекцию.
    :param data: данные из файла, которые были получены с серверов Pubchem.
    :param collection: Коллекция, в которую будут загружены данные.
    :return: Кортеж, где первый элемент – число загруженных документов, второй – незагруженных.
    """
    res = collection.insert_many(data, ordered=False).inserted_ids
    return len(res), (len(data) - len(res))


def delete_without_smiles(collection: pymongo.collection.Collection) -> int:
    """
    Для правильной работы базы необходимо удостовериться, что все документы в рабочей коллекции
    имеют поле ``CanonicalSMILES``. Документы не удовлетворяющие данному условию удаляются из
    коллекции.
    :param collection: Коллекция, в которой будет производиться операция удаления документов.
    :return: Количество удаленных документов.
    """
    count = collection.delete_many({"CanonicalSMILES": {"$exists": False}}).deleted_count
    return count


def sorting_search(smiles: str) -> List[str]:
    properties, molecules, mfp_counts = settings.get_collections
    q_mol = Chem.MolFromSmiles(smiles)
    substr_lst = substructure.SubSearch(q_mol, molecules, chirality=False)
    substr_set = set(substr_lst)
    res = []
    similar_search = similarity.SimSearchAggregate(q_mol, molecules, mfp_counts, 0.4)
    for score, smi in sorted(similar_search, reverse=True):
        res.append(smi)
        if smi in substr_set:
            substr_set.remove(smi)
    res.extend(substr_set)
    return res


def stages(result: List[str]):
    match = {
        "$match": {"rdkit_index": {"$in": result}}
    }
    group = {
        "$group": {
            "_id": 0,
            "avg_MolecularWeight": {"$avg": "$MolecularWeight"},
            "std_MolecularWeight": {"$stdDevPop": "$MolecularWeight"},
            "avg_XLogP": {"$avg": "$XLogP"},
            "std_XLogP": {"$stdDevPop": "$XLogP"},
            "avg_HBondDonorCount": {"$avg": "$HBondDonorCount"},
            "std_HBondDonorCount": {"$stdDevPop": "$HBondDonorCount"},
            "avg_HBondAcceptorCount": {"$avg": "$HBondAcceptorCount"},
            "std_HBondAcceptorCount": {"$stdDevPop": "$HBondAcceptorCount"},
            "avg_RotatableBondCount": {"$avg": "$RotatableBondCount"},
            "std_RotatableBondCount": {"$stdDevPop": "$RotatableBondCount"},
            "avg_AtomStereoCount": {"$avg": "$AtomStereoCount"},
            "std_AtomStereoCount": {"$stdDevPop": "$AtomStereoCount"},
            "avg_BondStereoCount": {"$avg": "$BondStereoCount"},
            "std_BondStereoCount": {"$stdDevPop": "$BondStereoCount"},
            "avg_Volume3D": {"$avg": "$Volume3D"},
            "std_Volume3D": {"$stdDevPop": "$Volume3D"},
        }
    }
    return match, group


def simple_search(smiles, skip, limit):
    properties, _, _ = settings.get_collections
    res = sorting_search(smiles)
    return properties.find({"rdkit_index": {"$in": res}}).skip(skip).limit(limit)


def summary_search(smiles):
    properties, _, _ = settings.get_collections
    res = sorting_search(smiles)
    pipeline = stages(res)
    res_summary = properties.aggregate(pipeline)
    res = simple_search(smiles, 0, 1)
    return res, res_summary


if __name__ == "__main__":
    smiles_ = "CCN1C=NC2=C(N=CN=C21)N"
    skip_ = 0
    limit_ = 10
    sorting_search(smiles_)

