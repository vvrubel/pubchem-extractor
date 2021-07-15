from typing import Any, Dict, Iterator, List, Tuple

import pymongo
from loguru import logger
from mongordkit.Database import registration
from pymongo.collection import Collection
from pymongo.database import Database
from rdkit import Chem

from molecad.errors import EmptySmilesError
from molecad.settings import settings
from molecad.validator import check_smiles


def prepare_db(db: Database, drop: bool) -> List[Collection]:
    """
    Если в базе данных есть коллекции, то она будет очищена, после чего коллекции будут созданы
    заново и на них будет создан индекс \"CID\".
    :param db: Объект базы данных.
    :param db: Флаг для очистки базы данных.
    :return: Список объектов коллекций.
    """
    lst = db.list_collection_names()
    if drop and len(lst) > 0:
        for item in lst:
            db.drop_collection(item)
            logger.info(f"Коллекция {item} удалена.", fg="yellow")

    names = [settings.pubchem, settings.molecules]
    collection_lst = []
    for name in names:
        collection: pymongo.collection.Collection = db.create_collection(name)
        collection.create_index("CID", unique=True)
        collection_lst.append(collection)
    return collection_lst


def register_from_smiles(smiles: str) -> Dict[str, Any]:
    """
    Создает объект схемы, определяющий модель документа молекул, и назначает индекс на коллекции.
    :param smiles: Строковое представление молекулы.
    :return: Схема документа молекулы.
    """
    rdmol = Chem.MolFromSmiles(smiles)
    scheme = registration.MolDocScheme()
    scheme.set_index("CanonicalSmiles")
    scheme.add_all_hashes()
    return scheme.generate_mol_doc(rdmol)


def populate_w_schema(data: List[Dict[str, Any]], mol_collection: Collection) -> int:
    n = 0
    for mol in data:
        smiles = mol["CanonicalSMILES"]
        cid = mol["CID"]
        try:
            check_smiles(smiles)
        except EmptySmilesError:
            continue
        else:
            scheme = register_from_smiles(smiles)
            scheme["CID"] = cid
            mol_collection.insert_one(scheme)
            n += 1
    return n


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


def retrieve_smiles(collection: pymongo.collection.Collection) -> Iterator[Dict[Any, Any]]:
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection: Коллекция, из которой будет извлекаться значение поля `CanonicalSMILES`.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    cursor = collection.find({}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
    for doc in cursor:
        yield doc.values()
