from typing import Any, Dict, List, Optional

import pymongo.errors
from loguru import logger
from pymongo import MongoClient

from molecad.settings import setup


class Connect:
    @staticmethod
    def get_connection():
        return MongoClient(
            host=setup.mongo_host,
            port=setup.mongo_port,
            username=setup.mongo_user,
            password=setup.mongo_password,
            authSource=setup.mongo_auth_source,
        )[setup.mongo_db_name]


db = Connect.get_connection()


def create_index(collection_name: str) -> None:
    """
    Создает на коллекции уникальный индекс по полю "CID".
    :param collection_name: Название коллекции, на которой будет создан индекс.
    :return: None.
    """
    try:
        db[collection_name].create_index([("CID", 1)], unique=True)
    except pymongo.errors.OperationFailure:
        logger.info(f"На коллекции {collection_name} индекс CID уже присутствует.")
    else:
        logger.info(f"На коллекции {collection_name} создан индекс – CID.")


def drop_collection(collection_name: str) -> None:
    """
    Если в коллекции есть документы, то она будет очищена.
    :param collection_name: Название коллекции, которая будет очищена.
    :return: None.
    """
    collection = db[collection_name]
    if collection.count_documents({}) > 0:
        collection.drop()
    logger.info(f"Коллекция {collection_name} очищена.")


def upload_data(data: List[Dict[str, Any]], collection_name: str) -> int:
    """
    Загружает данные в коллекцию и создает на ней индекс.
    :param data: взятые из файла или скачанные напрямую с серверов Pubchem.
    :param collection_name: Название коллекции, в которую будут загружены данные.
    :return: Число документов, загруженных за цикл работы функции.
    """
    # n = 0
    # for rec in data:
    #     try:
    #         db[collection_name].insert_one(rec)
    #     except pymongo.errors.DuplicateKeyError:
    #         pass
    #     else:
    #         n += 1
    # return n
    try:
        return db[collection_name].insert_many(data)
    except pymongo.errors.BulkWriteError:
        raise


def delete_without_smiles(collection_name: str) -> int:
    """
    В данных, загруженных с Pubchem, для некоторых молекул могут отсутствовать запрашиваемые
    поля - для наших задач критически важным является поле ``CanonicalSMILES``.
    Данная функция ищет в локальной базе документы, в которых это поле отсутствует, и удаляет их.
    :param collection_name: Название коллекции, в которой будет производиться операция удаления
    документов.
    :return: Количество удаленных документов.
    """
    collection = db[collection_name]
    count = collection.delete_many({"CanonicalSMILES": {"$exists": False}}).deleted_count
    logger.info(f"Удалено {count} документов")
    return count


def retrieve_smiles(collection_name: str, filters: Optional[Dict[str, Any]] = None):
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection_name: название коллекции, из которой будет извлекаться значение поля
    `CanonicalSMILES`.
    :param filters: Опциональный параметр, который может присутствовать, если необходимо
    установить поисковый фильтр.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    collection = db[collection_name]
    cursor = collection.find({}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
    smiles = {}
    for rec in cursor:
        smiles[rec["CID"]] = rec["CanonicalSMILES"]
    return smiles
