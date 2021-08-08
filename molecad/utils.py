import functools
import json
import time
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Tuple, TypeVar, Union

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
        print(
            f"Время выполнения {func.__name__!r}:\n "
            f"{hours} часов, {minutes} минут, {seconds:.2f} секунд."
        )
        return val

    return wrapper


def generate_ids(start: int, stop: int) -> Iterator[int]:
    for n in range(start, stop):
        yield n


def chunked(iterable: Iterable[T], maxsize: int) -> Iterable[List[T]]:
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield list(chunk)
            chunk.clear()
    if chunk:
        yield chunk


def concat(*args: Any, sep: str = "/") -> str:
    return sep.join(f"{i}" for i in args)


def parse_first_and_last(obj: List[Dict[str, Any]]) -> Tuple[int, int]:
    first = obj[0]["CID"]
    last = obj[-1]["CID"]
    return first, last


def create_dir(parent_dir: Path, first_id: int, last_id: int) -> Path:
    name = concat(first_id, last_id, sep="–")
    inc = 1
    new_dir = Path(parent_dir).resolve() / name
    while new_dir.exists():
        new_dir = Path(parent_dir).resolve() / concat(name, inc, sep="_")
        inc += 1

    new_dir.mkdir(parents=True, exist_ok=False)
    return new_dir


def file_path(dir_path: Path, first_id: int, last_id: int) -> Path:
    name = concat(first_id, last_id, sep="–")
    f_name = name + ".json"
    f_path = Path(dir_path) / f_name
    return f_path


def read_json(f_path: Path) -> Union[Dict[int, T], List[T]]:
    with open(f_path, "rt") as f:
        data = json.load(f)
        return data


def write_json(f_path: Path, data: Union[Dict[int, T], List[T]]) -> None:
    with open(f_path, "wt") as f:
        json.dump(data, f)


def converter(obj: Union[Dict[int, T], List[T]]) -> List[T]:
    if isinstance(obj, dict):
        return list(obj.values())
    else:
        return obj
