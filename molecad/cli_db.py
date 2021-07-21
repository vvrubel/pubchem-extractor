from typing import Any, Dict, List, Tuple

import click
from mongordkit.Database import registration
from pymongo.collection import Collection
from pymongo.database import Database
from rdkit import Chem


def drop_db(db: Database) -> None:
    """
    Если в базе данных есть коллекции, то она будет очищена, после чего коллекции будут созданы
    заново и на них будет создан индекс ``CID``.
    :param db: База данных, взятая из контекста ``click``.
    :return: None.
    """
    lst = db.list_collection_names()
    if len(lst) > 0:
        for item in lst:
            db.drop_collection(item)
            click.secho(f"Коллекция {item} удалена.", fg="red")


def create_indexes(*args: Collection) -> None:
    """
    Создает на переданных в функцию коллекциях уникальный индекс ``CID`` и индекс ``index``,
    а также индекс ``score`` на коллекции ``properties``.
    :param args: Список из коллекций ``properties`` и ``molecules``.
    :return: None.
    """
    for arg in args:
        arg.create_index("CID", unique=True)
        click.secho(f'На коллекции {arg.name} создан уникальный индекс "CID".', fg="green")
        arg.create_index("index")
        click.secho(f'На коллекции {arg.name} создан индекс – "index".', fg="green")
    args[0].create_index(["score", -1])
    click.secho(f'На коллекции {args[0].name} создан индекс – "score".', fg="green")


def register_from_smiles(smiles: str) -> Dict[str, Any]:
    """
    Создает объект схемы, определяющий модель документа молекул, и назначает индекс на коллекции.
    :param smiles: Строковое представление молекулы.
    :return: Словарь–репрезентация молекулы, сгенерированный из SMILES, который будет загружен в
    MongoDB в качестве документа в коллекцию ``molecules``.
    """
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise TypeError
    scheme = registration.MolDocScheme()
    scheme.set_index("CanonicalSmiles")
    return scheme.generate_mol_doc(rdmol)


def create_molecule(
    data: List[Dict[str, Any]], mol_collection: Collection
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Функция вызывается внутри команды ``populate`` и получает список из словарей (в будущем
    документов): в нем смотрит на наличие ключа ``CanonicalSMILES`` и, если он имеется, генерирует
    схему документа молекулы, добавляя в схему ключ ``CID``, вставляет документ в коллекцию
    ``molecules`` в случае успешного завершения предыдущих этапов. Из сгенерированного документа,
    достает значение по ключу ``index`` и добавляет его к ``data``. Также в ``data`` находит
    значение поля по ключу ``MolecularWeight`` и приводит его к типу ``float``, возвращая
    измененный словарь ``data`` из функции.
    :param data: Данные, загруженные из файла ``.json``.
    :param mol_collection: Коллекция, в которую будут вставлены документы, являющиеся
    представлением молекул.
    :return: Кортеж из измененной ``data`` и количества сгенерированных схем.
    """
    n = 0
    for molprop in data:
        if "CanonicalSMILES" not in molprop:
            continue
        smiles = molprop["CanonicalSMILES"]
        cid = molprop["CID"]
        try:
            scheme = register_from_smiles(smiles)
        except TypeError:
            # if scheme is None
            continue
        scheme["CID"] = cid

        mol_collection.insert_one(scheme)
        n += 1
        # из схемы берем поле "index" и добавляем его к данным
        molprop["index"] = scheme["index"]

        # приводим "MolecularWeight" к float
        mol_weight = molprop["MolecularWeight"]
        if mol_weight is not None:
            molprop["MolecularWeight"] = float(mol_weight)

    return data, n


def upload_data(data: List[Dict[str, Any]], collection: Collection) -> Tuple[int, int]:
    """
    Загружает данные в коллекцию.
    :param data: данные из файла, которые были получены с серверов Pubchem.
    :param collection: Коллекция, в которую будут загружены данные.
    :return: Кортеж, где первый элемент – число загруженных документов, второй – незагруженных.
    """
    res = collection.insert_many(data, ordered=False).inserted_ids
    return len(res), (len(data) - len(res))


def delete_broken(collection: Collection) -> Tuple[int, int]:
    """
    Для правильной работы базы необходимо удостовериться, что все документы в рабочей коллекции
    имеют поле ``CanonicalSMILES`` и сгенерированную схему. Документы не удовлетворяющие данным
    условиям удаляются из коллекции.
    :param collection: Коллекция, в которой будет производиться операция удаления документов.
    :return: Количество удаленных документов без SMILES и без схем.
    """
    without_smiles = collection.delete_many({"CanonicalSMILES": {"$exists": False}}).deleted_count
    without_schemes = collection.delete_many({"index": {"$exists": False}}).deleted_count
    return without_smiles, without_schemes
