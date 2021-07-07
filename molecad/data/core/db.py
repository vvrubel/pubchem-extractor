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
            authSource=setup.mongo_auth_source,
        )


client = Connect.get_connection()
db = client[setup.mongo_db_name]
collection = db[setup.mongo_db_collection]


def drop_collection() -> None:
    """
    Если в коллекции есть документы, то она будет очищена.
    """
    if collection.count_documents({}) > 0:
        collection.drop()
    logger.info(f"Коллекция {collection} очищена.")


def upload_data(data: List[Dict[str, Any]]) -> List[str]:
    """
    Загружает данные в коллекцию.
    :param data: взятые из файла или скачанные напрямую с серверов Pubchem.
    """
    return collection.insert_many(data).inserted_ids


def upload_test_data(data: List[Dict[str, Any]]) -> List[str]:
    """
    Загружает тестовые данные в тестовую коллекцию.
    :param data: взятые из файла или скачанные напрямую с серверов Pubchem.
    """
    test = db["test"]
    return test.insert_many(data).inserted_ids
