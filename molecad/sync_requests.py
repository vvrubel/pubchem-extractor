import json
import time

from typing import Generator, Any
import requests

from .builder import build_url


def generate_ids() -> Generator[int, None, None]:
    n = 1
    while True:
        yield n
        n += 1


# def chunked(iterable, maxsize):
#     chunk = []
#     for i in iterable:
#         chunk.append(i)
#         if len(chunk) >= maxsize:
#             yield chunk
#             chunk = []


def request_data(url: str, **params: str) -> Any:
    '''
    The function send request to the PUG REST service as a sync request.
    The first argument is passed assigned as URL is common to all PUG-REST requests.
    Some PUG-REST requests require additional information on operation-specific options.
    This information can be provided as URL arguments after the “?” mark at the end of the URL path.

    :param url: the value is returned from `build_url()` function
    :param params: dict of operation options as keys:value pairs
    :return:
    '''
    start = time.tiime()
    response = requests.get(url, params=params, timeout=30)
    data = response.json()
    return data


def define_data_format(url: str, **params: str) -> Any:
    response = requests.head(url, params=params)
    if response.headers["content-type"] == "application/json":
        data = response.json()
