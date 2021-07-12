from typing import Any, Dict, List, Optional

import pymongo.collection
import pymongo.errors
from loguru import logger
from mongordkit.Database import registration
from mongordkit.Search import substructure
from rdkit import Chem

from molecad.settings import setup

rdkit_schema = setup.get_db["test_rdkit_schema"]
mol = setup.get_db["molecules"]


def create_index(collection: pymongo.collection.Collection, index: str) -> None:
    """
    Создает на коллекции уникальный индекс по указанному полю.
    :param collection: Коллекция, на которой будет создан индекс.
    :param index: Строка, являющаяся полем в заданной коллекции, которое будет назначено в
    роли индекса.
    :return: None.
    """
    try:
        collection.create_index([(index, 1)], unique=True)
    except pymongo.errors.OperationFailure:
        logger.info(f"На коллекции {collection} индекс {index} уже присутствует.")
    else:
        logger.info(f"На коллекции {collection} создан индекс – {index}.")


def drop_collection(collection: pymongo.collection.Collection) -> None:
    """
    Если в коллекции есть документы, то она будет очищена.
    :param collection: Коллекция, которая будет очищена.
    :return: None.
    """
    if collection.count_documents({}) > 0:
        collection.drop()
    logger.info(f"Коллекция {collection} очищена.")


def upload_data(data: List[Dict[str, Any]], collection: pymongo.collection.Collection) -> int:
    """
    Загружает данные в коллекцию.
    :param data: данные из файла, которые были получены с серверов Pubchem.
    :param collection: Коллекция, в которую будут загружены данные.
    :return: Число документов, загруженных за цикл работы функции.
    """
    try:
        return collection.insert_many(data)
    except pymongo.errors.BulkWriteError:
        raise


def delete_without_smiles(collection: pymongo.collection.Collection) -> int:
    """
    Для правильной работы базы необходимо удостовериться, что все документы в рабочей коллекции
    имеют поле ``CanonicalSMILES``. Документы не удовлетворяющие данному условию удаляются из
    коллекции.
    :param collection: Коллекция, в которой будет производиться операция удаления документов.
    :return: Количество удаленных документов.
    """
    count = collection.delete_many({"CanonicalSMILES": {"$exists": False}}).deleted_count
    logger.info(f"Удалено {count} документов")
    return count


def retrieve_smiles(
    collection: pymongo.collection.Collection,
    filters: Optional[Dict[str, Any]] = None
) -> Dict[Any, Any]:
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection: Коллекция, из которой будет извлекаться значение поля `CanonicalSMILES`.
    :param filters: Опциональный параметр, который может присутствовать, если необходимо
    отфильтровать данные.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    cursor = collection.find({filters}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
    smiles = {}
    for rec in cursor:
        smiles[rec["CID"]] = rec["CanonicalSMILES"]
    return smiles


def schema_from_smiles(smiles: str) -> Dict[str, Any]:
    """
    Создает объект схемы, определяющий модель документа молекул, и назначает индекс на коллекции.
    :param smiles: Строковое представление молекулы.
    :return: Схема документа молекулы.
    """
    rdmol = Chem.MolFromSmiles(smiles)
    scheme = registration.MolDocScheme()
    scheme.set_index('CanonicalSmiles')
    return scheme.generate_mol_doc(rdmol)


def prepare_for_search():
    res = retrieve_smiles(mol)
    for cid, smiles in res.items():
        schema = schema_from_smiles(smiles)
        rdkit_schema.insert_one(schema)
        mol.update_one({"CID": cid}, {"$set": {"rdkit_schema_index": schema["index"]}}, upsert=True)
    create_index(rdkit_schema, "index")
    create_index(mol, "index")


def search(sub_smiles):
    q_mol = Chem.MolFromSmiles(sub_smiles)
    substructure.AddPatternFingerprints(rdkit_schema)
    index_lst = substructure.SubSearch(q_mol, rdkit_schema, chirality=False)
    # index_lst = substructure.SubSearchNaive(q_mol, rdkit_schema, chirality=False)
    cursor = mol.find({"rdkit_schema_index": {"$in": index_lst}})
    return list(cursor)


def main():
    prepare_for_search()
    res = search("CC(=O)OC")
    print([i['CanonicalSMILES'] for i in res])


if __name__ == "__main__":
    main()
