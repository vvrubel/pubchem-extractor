from typing import Any, Dict, List, Optional, Tuple

import pymongo
from loguru import logger
from mongordkit.Database import registration
from mongordkit.Search import substructure
from rdkit import Chem

from molecad.errors import CreateIndexError
from molecad.validator import is_index_exist
from molecad.settings import setup


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
        logger.warning(f"На коллекции {collection.name} индекс {index} уже создан.")
        raise CreateIndexError
    else:
        collection.create_index(index, unique=True)
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
    logger.info(f"Удалено {count} документов")
    return count


def retrieve_smiles(collection: pymongo.collection.Collection) -> Dict[Any, Any]:
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection: Коллекция, из которой будет извлекаться значение поля `CanonicalSMILES`.
    :param filters: Опциональный параметр, который может присутствовать, если необходимо
    отфильтровать данные.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    cursor = collection.find({}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
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


def prepare_for_search(
    collection: pymongo.collection.Collection,
    rdkit_schema: pymongo.collection.Collection
) -> None:
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
