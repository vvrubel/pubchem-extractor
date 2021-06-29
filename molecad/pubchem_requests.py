import requests
import time
from loguru import logger
from typing import (
    Any,
    Dict,
    Generator,
    Iterable,
    TypeVar
)

from .url_builder import prepare_request
from .types_ import (
    Domain,
    NamespCmpdAlone,
    OperationAlone,
    OperationComplex,
    Out,
    PropertyTags,
)
from .data_handler import save_data_json

T = TypeVar("T")


def generate_ids(
        start: int = 1,
        stop: int = 201
) -> Generator[int, None, None]:
    """
    Simple generator of CID values.
    :param start: default value is 1, but for custom downloads can take any
    value - usually is ``stop`` from the previous download.
    :param stop: any value below 156 millions.
    :return: generator of integers.
    """
    for n in range(start, stop):
        yield n


def chunked(
        iterable: Iterable[T],
        maxsize: int
) -> Generator[list[T], None, None]:
    """
    Takes an iterable with definite type and divides it to chunks with equal
    size.
    :param iterable: sequence passed from ``generate_ids()`` function call.
    :param maxsize: max number of elements in a chunk.
    :return: yielded chunks must be passed to ``join_w_comma()`` before it
    will be put into ``input_specification()`` function.
    """
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
    """
    PubChem REST service has the limitations:
    No more than five requests per second.
    No more than 400 requests per minute.
    No more than 300 second running time on PubChem servers per minute.
    :param iterable: sequence of ids or chunks.
    :param waiting_time: 60 sec.
    :param maxsize: 400 requests.
    :return: generates items according with time limit of requests.
    """
    window = []
    for i in iterable:
        yield i
        t = time.monotonic()
        window.append(t)
        while t - waiting_time > window[0]:
            window.pop(0)
        if len(window) > maxsize:
            t0 = window[0]
            delay = t - t0
            time.sleep(delay)


def request_data_json(url: str, **params: str) -> Dict[str, Any]:
    """
    Function sends a synchronous "GET" request to the PUG REST service.
    The first argument - ``url`` is common to all requests; the additional
    arguments are required for some operation-specific options.
    If you were constructing the URL by hand, this data would be given as
    ``key:value`` pairs in the URL after ``?`` mark at the end.
    :param url: returned from ``build_url()`` function.
    :param params: operation options as a dict of keys and values which are
    string type.
    :return: response in JSON format.
    """
    response = requests.get(url, params=params).json()
    return response


def main():
    domain = Domain.COMPOUND
    namespace = NamespCmpdAlone.CID
    operation = OperationComplex.PROPERTY
    tags = (
        PropertyTags.MOLECULAR_FORMULA,
        PropertyTags.MOLECULAR_WEIGHT,
        PropertyTags.CANONICAL_SMILES,
        PropertyTags.INCHI,
        PropertyTags.IUPAC_NAME,
        PropertyTags.XLOGP,
        PropertyTags.H_BOND_DONOR_COUNT,
        PropertyTags.H_BOND_ACCEPTOR_COUNT,
        PropertyTags.ROTATABLE_BOND_COUNT,
        PropertyTags.ATOM_STEREO_COUNT,
        PropertyTags.BOND_STEREO_COUNT,
        PropertyTags.VOLUME_3D,
    )
    output = Out.JSON

    results = {}
    t_start = time.monotonic()
    chunks = chunked(generate_ids(), 100)
    for i in delay_iterations(chunks, 60.0, 400):
        url = prepare_request(domain, namespace, i, operation, output, tags)
        logger.debug("Requesting URL: {}", url)
        try:
            res = request_data_json(url)
        except requests.HTTPError:
            logger.error("Error occurred: {}", exc_info=True)
            break
        else:
            logger.debug("Response content: {}", res)
            for k, v in zip(i, res['PropertyTable']['Properties']):
                results[k] = v
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    logger.info("Download took {}", t_run)

    res_file = save_data_json(results)
    logger.info("Data saved in {}", res_file)


if __name__ == '__main__':
    main()
