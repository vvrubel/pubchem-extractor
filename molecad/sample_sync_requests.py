import requests
import time
from typing import (
    Any,
    Dict,
    Generator,
    Iterable,
    List
)

from molecad.types_ import T


def generate_ids(
        start: int = 1,
        stop: int = 500001
) -> Generator[int, None, None]:
    for n in range(start, stop):
        yield n


def chunked(
        iterable: Iterable[T],
        maxsize: int
) -> Generator[List[T], None, None]:
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield chunk
            chunk = []


def delay_iterations(
        iterable: Iterable[T],
        waiting_time: float,
        maxsize: int
) -> Generator[T, None, None]:
    window = []
    for i in iterable:
        yield i
        t = time.monotonic()
        window.append(t)
        while t - waiting_time > window[0]:
            del window[0]
        if len(window) > maxsize:
            t0 = window[0]
            delay = t - t0
            time.sleep(delay)


def request_data_json(url: str, **params: str) -> Dict[str, Any]:
    response = requests.get(url, params=params)
    data = response.json()
    return data
