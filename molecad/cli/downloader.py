import time
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    TypeVar,
    Union,
)

import requests
from loguru import logger

from molecad.cli.errors import BadDomainError, BadNamespaceError, BadOperationError
from molecad.cli.models import Domain, NamespCmpd, OperationComplex, Out, PropertyTags
from molecad.cli.validator import (
    is_complex_operation,
    is_compound,
    is_namespace_search,
    is_simple_namespace,
    is_simple_operation,
)
from molecad.utils import chunked, concat, generate_ids, parse_first_and_last

IdT = TypeVar("IdT", int, str)
T = TypeVar("T")
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def url_builder(
    ids: Union[IdT, Iterable[IdT]],
    *,
    domain: str,
    namespace_prefix: str,
    namespace_suffix: Optional[str] = None,
    operation: str,
    tags: Optional[Sequence[str]] = None,
    output: str,
) -> str:
    if not is_compound(domain):
        raise BadDomainError

    if is_simple_namespace(namespace_prefix, namespace_suffix):
        joined_namespaces = namespace_prefix
    elif is_namespace_search(namespace_prefix, namespace_suffix):
        joined_namespaces = concat(namespace_prefix, namespace_suffix)
    else:
        raise BadNamespaceError

    joined_identifiers = concat(*ids, sep=",")
    input_spec = concat(domain, joined_namespaces, joined_identifiers)

    if is_complex_operation(operation, tags):
        joined_tags = concat(*tags, sep=",")
        operation_spec = concat(operation, joined_tags)
    elif is_simple_operation(operation, tags):
        operation_spec = operation
    else:
        raise BadOperationError

    url = concat(BASE_URL, input_spec, operation_spec, output)
    return url


def request_data_json(url: str, **operation_options: str) -> List[Dict[str, Any]]:
    response = requests.get(url, params=operation_options)
    response.raise_for_status()
    res = response.json()
    return res["PropertyTable"]["Properties"]


def delay_iterations(
    iterable: Iterable[T], waiting_time: float = 60.0, maxsize: int = 400
) -> Iterator[T]:
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


def execute_requests(start: int, stop: int, maxsize: int = 100) -> Iterator[Dict[str, Any]]:
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
    d = {
        "domain": Domain.COMPOUND,
        "namespace_prefix": NamespCmpd.CID,
        "namespace_suffix": None,
        "operation": OperationComplex.PROPERTY,
        "tags": tags,
        "output": Out.JSON,
    }
    chunks = chunked(generate_ids(start, stop), maxsize)
    for chunk in delay_iterations(chunks):
        try:
            url = url_builder(chunk, **d)
            res = request_data_json(url)
        except requests.exceptions.HTTPError:
            logger.opt(exception=True).error("Ошибка при выполнении запроса.")
            break
        except KeyError:
            logger.opt(exception=True).error("Неверный формат ответа от сервера.")
            break
        except BadDomainError:
            logger.error("В данной версии сервиса поиск возможен только по базе данных 'Compound'.")
            break
        except BadNamespaceError:
            logger.error("Пространство имен поиска по базе данных 'Compound' задано некорректно.")
            break
        except BadOperationError:
            logger.error("Ошибка при составлении операции.")
            break
        else:
            first, last = parse_first_and_last(res)
            logger.info(f"Получены CID: {first}–{last}")
            for rec in res:
                yield rec
