from typing import Any, Dict, List, Tuple

import pymongo

from .errors import CreateIndexError
from .validator import is_index_exist


def create_index(collection: pymongo.collection.Collection, index: str) -> None:
    """
    Пробует создать на коллекции уникальный индекс по указанному полю и говорит пользователю,
    что он создан в случае успеха; если такой индекс уже задан на коллекции, говорит, что индекс
    создать не получилось, так как он уже присутствует.
    :param collection: Коллекция, на которой будет создан индекс.
    :param index: Поле-строка, которое будет назначено в качестве индекса.
    :return: None.
    """
    if is_index_exist(collection, index):
        raise CreateIndexError
    else:
        collection.create_index(index, unique=True)


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


def retrieve_smiles(collection: pymongo.collection.Collection) -> Dict[Any, Any]:
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection: Коллекция, из которой будет извлекаться значение поля `CanonicalSMILES`.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    cursor = collection.find({}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
    smiles = {}
    for rec in cursor:
        smiles[rec["CID"]] = rec["CanonicalSMILES"]
    return smiles
