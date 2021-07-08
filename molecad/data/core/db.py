from typing import Any, Dict, List

import pymongo.results
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


def use_collection(collection_name: str):
    """
    :param collection_name: название запрашиваемой коллекции.
    :return: Возвращает коллекцию по запрашиваемому имени.
    """
    return db[collection_name]


def drop_collection(collection_name: str) -> None:
    """
    Если в коллекции есть документы, то она будет очищена.
    """
    collection = use_collection(collection_name)
    if collection.count_documents({}) > 0:
        collection.drop()
    logger.info(f"Коллекция {collection} очищена.")


def upload_data(data: List[Dict[str, Any]], collection_name: str) -> None:
    """
    Загружает данные в коллекцию.
    :param data: взятые из файла или скачанные напрямую с серверов Pubchem.
    :param collection_name: название коллекции, в которую будут загружены данные.
    """
    collection = use_collection(collection_name)
    return collection.insert_many(data)


def retrieve_smiles(collection_name: str, filters):
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов
    :param collection_name: название коллекции, из которой будет извлекаться smiles.
    :param filters:
    :return:
    """
    collection = use_collection(collection_name)
    cursor = collection.find(filters, {"CanonicalSMILES": 1, "_id": 0})
    smiles = []
    for n in cursor:
        smiles.append(n)
    return smiles


def delete_items_without_smiles(collection_name: str):
    """
    В данных, загруженных с Pubchem, для некоторых молекул могут отсутствовать некоторые из
    запрашиваемых полей - для наших задач критически важным является поле ``CanonicalSMILES``.
    Данная функция ищет в локальной базе документ, в которых это поле отсутствует, и удаляет их.
    :param collection_name: название коллекции, в которой будет производиться операция удаления
    документов.
    :return: количество удаленных документов.
    """
    collection = use_collection(collection_name)
    collection.delete_many({"CanonicalSMILES": {"$exists": "false"}})
    cnt = pymongo.results.DeleteResult.deleted_count
    logger.info(f"{cnt} document were deleted")
