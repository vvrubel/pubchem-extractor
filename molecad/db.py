from typing import Any, Dict, Iterator, List, Tuple

import pymongo
import rdkit.Chem
from loguru import logger
from mongordkit.Database import registration
from mongordkit.Search import similarity, substructure
from pymongo.collection import Collection
from pymongo.database import Database
from rdkit import Chem

from molecad.errors import EmptySmilesError
from molecad.settings import settings
from molecad.validator import check_smiles

test_smiles = "S(=O)(=O)NC(=O)N"


def prepare_db(db: Database, drop: bool) -> List[Collection]:
    """
    Если в базе данных есть коллекции, то она будет очищена, после чего коллекции будут созданы
    заново и на них будет создан индекс \"CID\".
    :param db: Объект базы данных.
    :param drop: Флаг для очистки базы данных.
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


def populate_w_schema(
    data: List[Dict[str, Any]], mol_collection: Collection
) -> Tuple[List[Dict[str, Any]], int]:
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


def retrieve_smiles(collection: pymongo.collection.Collection) -> Iterator[Dict[Any, Any]]:
    """
    Функция позволяет извлечь поле `CanonicalSMILES` из документов.
    :param collection: Коллекция, из которой будет извлекаться значение поля `CanonicalSMILES`.
    :return: Строка, являющаяся представлением молекулы в формате SMILES.
    """
    cursor = collection.find({}, {"CanonicalSMILES": 1, "CID": 1, "_id": 0})
    for doc in cursor:
        yield doc.values()


def substructure_search(
    smiles: str,
    molecules: pymongo.collection.Collection,
) -> Tuple[rdkit.Chem.Mol, List[str]]:
    q_mol = Chem.MolFromSmiles(smiles)
    index_lst = substructure.SubSearch(q_mol, molecules, chirality=False)
    return q_mol, index_lst


def similarity_search(smiles, db, molecules, mfp_counts, permutations):
    q_mol = Chem.MolFromSmiles(smiles)
    res1 = similarity.SimSearch(q_mol, molecules, mfp_counts, 0.4)
    print(f"similaritySearch: {res1}")
    res2 = similarity.SimSearchAggregate(q_mol, molecules, mfp_counts, 0.4)
    print(f"similaritySearchAggregate: {res2}")
    res3 = similarity.SimSearchLSH(q_mol, db, molecules, permutations, mfp_counts, threshold=0.8)
    print(f"similaritySearchLSH: {res3}")
    results = [res1, res2, res3]
    return results


# TODO in next commit
def search(smiles, skipped, limited):
    db = settings.get_db
    q_mol, index_lst = substructure_search(smiles, settings.molecules)
    cursor = settings.pubchem.find({"rdkit_index": {"$in": index_lst}})
    return cursor.skip(skipped).limit(limited)


# TODO in next commit
if __name__ == "__main__":
    collection, mol_collection, mfp_collection, perm_collections = settings.get_collections
    substructure_search(test_smiles, mol_collection)
    # main(test_smiles, 0, 5)
