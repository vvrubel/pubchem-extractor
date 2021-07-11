import json
import uuid
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, TypeVar, Union

from loguru import logger

from molecad.errors import DirExistsError

T = TypeVar("T")


def generate_ids(start: int, stop: int) -> Iterator[int]:
    """
    Простой генератор значений CID.
    :param start: любое целое положительное число такие, что ``start < stop``.
    :param stop: любое целое положительное число такие, что ``start < stop``.
    :return: генератор целых положительных чисел.
    """
    for n in range(start, stop):
        yield n


def chunked(iterable: Iterable[T], maxsize: int) -> Iterable[List[T]]:
    """
    Принимает итерируемый объект со значениями одинакового типа и делит его на списки одинаковой
    длины, равной ``maxsize``, последний чанк содержит оставшиеся элементы.
    :param iterable: итерируемый объект.
    :param maxsize: максимальный размер чанка.
    :return: выбрасывает чанки заданной длины.
    """
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield list(chunk)
            chunk.clear()
    if chunk:
        yield chunk


def join_w_comma(*args: T) -> str:
    """
    Функция принимает на вход последовательность аргументов, приводит каждый из них к строке,
    после чего конкатенирует их с помощью запятой.
    :param args: любая последовательность аргументов одинакового типа.
    :return: строка разделенных запятой и без пробела значений.
    """
    return ",".join(f"{i}" for i in args)


def concat(*args: T, sep="/") -> str:
    """
    Функция принимает на вход последовательность аргументов, приводит каждый из них к строке,
    после чего конкатенирует их с помощью строки, переданной в `sep`.
    :param args: последовательность аргументов одинакового типа.
    :param sep: строка, являющаяся разделителем, с помощью которой будут соединены аргументы,
    переданные в `args`. Если значение не определено, то по умолчанию будет использоваться "/".
    :return: строка, соединенная с помощью `sep`.
    """
    return sep.join(f"{i}" for i in args)


def check_dir(dir_path: Path, start_id: int, stop_id: int) -> Path:
    """
    Рекурсивная функция, которая пробует создать поддиректорию c именем `name` в директории
    `dir_path`, если такая поддиректория уже существует, то выбрасывает ошибку.
    :param dir_path: Путь до поддиректории.
    :param start_id: Первое значение из запрашиваемых идентификаторов.
    :param stop_id: Последнее значение из запрашиваемых идентификаторов.
    :return: Путь до созданной директории.
    """

    name = concat(start_id, stop_id, sep="–")
    try:
        new_dir = Path(dir_path).resolve() / str(name)
        new_dir.mkdir(parents=True, exist_ok=False)
    except DirExistsError:
        logger.info("Директория уже существует.")
        raise
    else:
        return new_dir


def file_name(dir_path: Path) -> Path:
    """
    Придумывает имя для файла в формате json.
    :param dir_path: имя папки в которую будет в последующем записан файл.
    :return: путь до файла.
    """
    name = str(uuid.uuid4()) + ".json"
    f_path = Path(dir_path) / name
    return f_path


def read(f_path: Path) -> Union[Dict[int, T], List[T]]:
    """
    Читает данные из файла.
    :param f_path: Абсолютный путь до файла.
    :return: JSON объект.
    """
    with open(f_path, "rt") as f:
        data = json.load(f)
        return data


def write(
    f_path: Path,
    data: Union[Dict[int, T], List[T]],
) -> None:
    """
    Пишет данные в файл.
    :param f_path: Абсолютный путь до файла.
    :param data: JSON объект.
    :return: None.
    """
    with open(f_path, "wt") as f:
        json.dump(data, f)


def converter(obj: Union[Dict[int, T], List[T]]) -> List[T]:
    """
    Согласует тип данных объекта.
    :param obj: может иметь тип словаря или списка.
    :return: в случае если объект имеет тип словаря, то функция возвращает список его значений;
    если же объект имеет тип списка (оставшиеся случаи), то функция возвращает сам объект.
    """
    if isinstance(obj, dict):
        return list(obj.values())
    else:
        return obj
