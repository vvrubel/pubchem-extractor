import json
import uuid
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, TypeVar, Union

from loguru import logger

from molecad.errors import DirExistsError

T = TypeVar("T")


def generate_ids(start: int = 1, stop: int = 201) -> Iterator[int]:
    """
    Простой генератор значений CID.
    .. note:: В рамках формирования базы данных интервал идентификаторов был равен (1, 500001).
    :param start: по умолчанию равно 1, но может быть заменено на любое положительное число -
    для скачивания порциями равно значению ``stop`` в предыдущей порции загрузки.
    :param stop: любое значение до 156 миллионов; в тестовом режиме установлено значение 201 для
    получения быстрого результата.
    :return: генератор целых положительных чисел.
    """
    for n in range(start, stop):
        yield n


def chunked(iterable: Iterable[T], maxsize: int) -> Iterable[List[T]]:
    """
    Принимает итерируемый объект со значениями одинакового типа и делит его на чанки одинаковой
    длины, равной ``maxsize``.
    :param iterable: последовательность, пришедшая из функции ``generate_ids``.
    :param maxsize: максимальное число элементов в чанке, для формирования базы данных равно
    1000, тестовые запросы должны выполняться со значением 100.
    :return: выбрасывает списки элементов - чанки, которые затем необходимо передать в функцию
    ``join_w_comma()``, после чего полученная строка может быть передана в ``input_specification``.
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
    Функция форматирует входящую последовательность аргументов, соединяя элементы запятой без
    пробелов и приводит полученное к строке.
    :param args: любая последовательность аргументов одинакового типа.
    :return: строка разделенных запятой и без пробела значений.
    """
    return ",".join(f"{i}" for i in args)


def check_dir(dir_path: Path) -> None:
    """
    Пробует создать директорию, если директория существует по указанному пути,
    то кидает ошибку и просит указать другое имя.
    :param dir_path: путь до несуществующей директории.
    """
    try:
        dir_path.mkdir(parents=True, exist_ok=False)
    except DirExistsError:
        logger.info("Директория уже существует. Создайте поддиректорию.")
        raise


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
