import functools
import json
import time
import urllib.parse
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Tuple, TypeVar, Union

from loguru import logger

T = TypeVar("T")


def timer(func: Callable[..., T]) -> Callable[..., T]:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.monotonic()
        val = func(*args, **kwargs)
        end = time.monotonic()
        work_time = end - start
        hours = work_time // 3600
        minutes = (work_time % 3600) // 60
        seconds = (work_time % 3600) % 60
        logger.success(
            f"Время выполнения {func.__name__!r}:\n "
            f"{hours} часов, {minutes} минут, {seconds:.2f} секунд."
        )
        return val

    return wrapper


def generate_ids(start: int, stop: int) -> Iterator[int]:
    """
    Простой генератор значений 'CID'.
    :param start: Любое целое положительное число такие, что ``start < stop``.
    :param stop: Любое целое положительное число такие, что ``start < stop``.
    :return: Генератор целых положительных чисел.
    """
    for n in range(start, stop):
        yield n


def chunked(iterable: Iterable[T], maxsize: int) -> Iterable[List[T]]:
    """
    Принимает итерируемый объект со значениями одинакового типа и делит его на контейнеры (чанки)
    одинаковой длины, равной ``maxsize``, последний контейнер содержит оставшиеся элементы.
    :param iterable: Итерируемый объект.
    :param maxsize: Максимальный размер контейнера.
    :return: Выбрасывает чанки заданной длины.
    """
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield list(chunk)
            chunk.clear()
    if chunk:
        yield chunk


def concat(*args: Any, sep="/") -> str:
    """
    Функция принимает на вход последовательность аргументов, приводит каждый из них к строке,
    после чего конкатенирует их с помощью строки, переданной в ``sep``.
    :param args: Последовательность аргументов одинакового типа.
    :param sep: Строка, являющаяся разделителем, с помощью которой будут соединены аргументы,
    переданные в ``args``. Если значение не определено, то по умолчанию будет использоваться ``/``.
    :return: Строка, соединенная из строковых представлений элементов ``args`` с помощью ``sep``.
    """
    return sep.join(f"{i}" for i in args)


def parse_first_and_last(obj: List[Dict[str, Any]]) -> Tuple[int, int]:
    """
    Из первого и последнего элементов списка, представленного последовательностью словарей,
    извлекает значение поля словаря по ключу ``CID``.
    :param obj: Список из словарей, которые имеют ключ ``CID``.
    :return: Кортеж со значениями поля ``CID`` для первого и последнего элементов списка.
    """
    first = obj[0]["CID"]
    last = obj[-1]["CID"]
    return first, last


def create_dir(parent_dir: Path, first_id: int, last_id: int) -> Path:
    """
    Рекурсивная функция, которая пробует создать поддиректорию c именем ``name`` в директории
    ``dir_path``, если такая поддиректория уже существует, то прибавляет к ``name`` суффикс и
    переходит в начало цикла проверки. При отсутствии директории с заданным именем – создает её.
    Имя поддиректории генерируется из значений запрашиваемых идентификаторов.
    :param parent_dir: Путь до родительской директории.
    :param first_id: Первое значение из запрашиваемых идентификаторов.
    :param last_id: Последнее значение из запрашиваемых идентификаторов.
    :return: Путь до созданной поддиректории.
    """
    name = concat(first_id, last_id, sep="–")
    inc = 1
    new_dir = Path(parent_dir).resolve() / name
    while new_dir.exists():
        new_dir = Path(parent_dir).resolve() / concat(name, inc, sep="_")
        inc += 1

    new_dir.mkdir(parents=True, exist_ok=False)
    return new_dir


def file_path(dir_path: Path, first_id: int, last_id: int) -> Path:
    """
    Генерирует имя файла в формате ``.json``.
    :param dir_path: Путь до директории, в которую будет сохранен файл.
    :param first_id: Первое значение из сохраняемых идентификаторов.
    :param last_id: Последнее значение из сохраняемых идентификаторов.
    :return: Путь до файла.
    """
    name = concat(first_id, last_id, sep="–")
    f_name = name + ".json"
    f_path = Path(dir_path) / f_name
    return f_path


def read_json(f_path: Path) -> Union[Dict[int, T], List[T]]:
    """
    Читает данные из файла.
    :param f_path: Абсолютный путь до файла.
    :return: JSON объект.
    """
    with open(f_path, "rt") as f:
        data = json.load(f)
        return data


def write_json(f_path: Path, data: Union[Dict[int, T], List[T]]) -> None:
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
    :param obj: Может иметь тип словаря или списка.
    :return: В случае если объект имеет тип словаря, то функция возвращает список его значений;
    если же объект имеет тип списка (оставшиеся случаи), то функция возвращает сам объект.
    """
    if isinstance(obj, dict):
        return list(obj.values())
    else:
        return obj


def url_encoder(route: str, query: Dict[str, Union[str, int]]) -> str:
    """
    Энкодер URL-адреса для api. Нужен для того, чтобы не использовать сторонние сервисы для
    генерирования URL.
    :param route: Путь до ручки api, например: '/v1/compound'.
    :param query: Параметры, передаваемые в URL-адрес в формате словаря.
    Например: ``query = {'smiles': 'NC(=O)N', 'skip': 0, 'limit': 10}``.
    :return: Строка URL-адреса.
    """
    params = urllib.parse.urlencode(query)
    return f"http://127.0.0.1:8000{route}?{params}"
