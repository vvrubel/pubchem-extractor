from typing import Any, Dict, List

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
            authSource=setup.mongo_db_name,
        )


client = Connect.get_connection()
db = client[setup.mongo_db_name]


def use_collection(collection_name):
    """
    :param collection_name: название запрашиваемой коллекции.
    :return: Возвращает коллекцию по запрашиваемому имени.
    """
    return db[collection_name]


def drop_collection(name) -> None:
    """
    Если в коллекции есть документы, то она будет очищена.
    """
    collection = use_collection(name)
    if collection.count_documents({}) > 0:
        collection.drop()
    logger.info(f"Коллекция {collection} очищена.")


def upload_data(data: List[Dict[str, Any]], name: str) -> List[str]:
    """
    Загружает данные в коллекцию.
    :param data: взятые из файла или скачанные напрямую с серверов Pubchem.
    :param name: название коллекции, в которую будут загружены данные.
    """
    collection = use_collection(name)
    return collection.insert_many(data)
