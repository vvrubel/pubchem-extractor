import requests
import time
from typing import (
    Any,
    Dict,
    Generator,
    Iterable,
    List,
)


def generate_ids(
        start: int = 1,
        stop: int = 501
) -> Generator[int, None, None]:
    """
    Simple generator CID values.
    :param start: default value is 1, but for the next downloads may be
    set up to any value from the last call value - usually from the ``stop`` .
    :param stop: any value below 156 millions.
    :return: iterable sequence of integers.
    """
    for n in range(start, stop):
        yield n


def chunked(
        iterable: Iterable[int],
        maxsize: int
) -> Generator[List[int], None, None]:
    """
    Takes sequence and divides it to chunks with equal size.
    :param iterable: inputted sequence from ``generate_ids()`` function call.
    :param maxsize: size of chunk.
    :return: yielded chunks must be passed to ``join_w_comma()`` before it
    will be put into ``input_specification()``.
    """
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield chunk
            chunk = []


def delay_iterations(
        iterable: Iterable[int],
        waiting_time: float,
        maxsize: int
) -> Generator[int, None, None]:
    """
    PubChem REST service has the limitations:
    No more than five requests per second.
    No more than 400 requests per minute.
    No more than 300 second running time on PubChem servers per minute.
    :param iterable: ids or chunks of ids.
    :param waiting_time: 60 sec.
    :param maxsize: 400 requests.
    :return: ids or chunks according with time limit of requests.
    """
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
    """
    he function sends a synchronous "GET" request to the PUG REST service.
    The first argument - ``url``, is common to requests; the additional
    arguments are required for some operation-specific options.
    If you were constructing the URL by hand, this data would be given as
    ``key:value`` pairs in the URL after a ``?`` mark at the end.
    :param url: string value was returned from ``build_url()`` function.
    :param params: operation specification as a dict with string keys and
    values.
    :return: response in JSON format.
    """
    response = requests.get(url, params=params)
    data = response.json()
    return data
