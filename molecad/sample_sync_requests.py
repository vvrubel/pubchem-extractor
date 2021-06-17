import requests
import time
from typing import (
    Any,
    Generator,
    Iterable,
    List,
    TypeVar
)

T = TypeVar("T")


def generate_ids() -> Generator[int, None, None]:
    """
    Generates CID - the parameter of identifier in input specification
    :return: CID's generator
    """
    n = 1
    while True:
        yield n
        n += 1


def chunked(iterable: Iterable[T], maxsize: int) -> Generator[List[T], None, None]:
    """
    For faster requests it groups the CIDs into containers with definite size
    :param iterable: generated CIDs
    :param maxsize: size og the container
    :return: containers with CIDs
    """
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield chunk
            chunk = []


def delay_iterations(iterable: Iterable[T], width: float, maxsize: int) -> Generator[T, None, None]:
    """
    Limitations:
    No more than five requests per second.
    No more than 400 requests per minute.
    No more than 300 second running time on PubChem servers per minute.
    :param iterable: chunks or ids
    :param width: time window for counting requests
    :param maxsize: limit number of requests
    :return: chuck or id and sleeps if limit of requests is over
    """
    window = []
    for i in iterable:
        yield i
        t = time.monotonic()
        window.append(t)
        while t - width > window[0]:
            del window[0]
        if len(window) > maxsize:
            t0 = window[0]
            delay = t - t0
            time.sleep(delay)


def request_data_json(url: str, **params: str) -> Any:
    """
    The function sends a synchronous request to the PUG REST service.
    The first argument is common to all PUG-REST requests, others require operation-specific options.
    The last can be provided as a dict using the `params` keyword arguments.
    It takes the same result as you construct the URL by hands putting key/value pairs after the “?” mark
    at the end of the URL path.
    :param url: string value is returned from `build_url()` function
    :param params: dict of operation options as keys:value pairs
    :return: response in JSON format
    """
    response = requests.get(url, params=params)
    data = response.json()
    return data


def execute_request(url: str, params: dict[str, Any]) -> Any:
    print(url)
    res = request_data_json(url, **params)
    return res


# TODO
def define_data_format(url: str, **params: str) -> Any:
    response = requests.head(url, params=params)
    if response.headers["content-type"] == "application/json":
        data = response.json()
        return data
    if response.headers["content-type"] == "image/png":
        with open("new_image", "rb") as f:
            pass
