import time
from typing import (
    Any, Dict, Iterable, Iterator, List, Optional, Sequence, TypeVar, Union,
)

import requests
from loguru import logger

from .downloader_types import Domain, NamespCmpd, OperationComplex, Out, PropertyTags
from .errors import BadDomainError, BadNamespaceError, BadOperationError
from .utils import chunked, concat, generate_ids
from .validator import (
    is_complex_operation,
    is_compound,
    is_namespace_search,
    is_simple_namespace,
    is_simple_operation,
)

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
    """
    Команда генерирует строку URL-адреса, необходимую для выполнения запроса в службу PubChem.
    :param ids: Последовательность из значений одинакового типа, которые интерпретируются как
    идентификаторы соединений; последовательность идентификаторов может быть получена в
    результате выполнения функции ``chunked``.
    :param domain: Принимает значение из класса ``Domain``. В текущей версии сервиса доступен
    поиск только по базе данных ``Compound``, значение по умолчанию = ``Domain.COMPOUND``.
    :param namespace_prefix: Является первой частью, описывающей пространство имен
    поиска и может принимать значения из классов ``NamespCmpd`` и ``PrefixSearch``. Значение по
    умолчанию = ``NamespCmpd.CID``.
    :param namespace_suffix: Является второй частью, описывающей пространство имен
    поиска и может принимать значения из класса ``SuffixSearch`` или None; по умолчанию = None.
    :param operation: Является частью URL-адреса, которая описывает, какие действия необходимо
    выполнить над запрашиваемыми идентификаторами. Может принимать значения из классов ``Operation``
    и ``OperationComplex``. По умолчанию = ``OperationComplex.PROPERTY``.
    :param tags: Строка содержит последовательность из тегов, принадлежащих классу
    ``PropertyTags``; определяется только в случае если ``operation`` принадлежит к классу
    ``OperationComplex.PROPERTY``, иначе равно None.
    :param output: Принимает значение из класса ``Out``; по умолчанию - ``Out.JSON``.
    :return: Если все параметры принадлежат к указанным классам, то генерируется URL-адрес, который
    используется для запроса в базу данных Pubchem; иначе кидается соответствующая ошибка.
    """

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
    """
    Функция добавляет параметры к URL, если они переданы в переменную ``operation_options`` и
    отправляет синхронный GET-запрос к серверам PUG REST баз данных Pubchem.
    Если бы вы создавали URL-адрес вручную, пары значений из словаря ``operation_options``
    записывались бы в конец URL-адреса после знака ``?`` в виде ``key=value``, а при наличии
    более одной пары последние объединялись знаком ``&``.
    :param url: Возвращается из функции ``url_builder`` и является обязательным для всех
    запросов.
    :param operation_options: Может быть None или словарем из строк в случае определенных
    значений параметра ``operation`` в функции ``url_builder``.
    :return: Ответ от сервера в формате JSON, содержимое которого является списком из словарей.
    """
    try:
        response = requests.get(url, params=operation_options).json()
        return response["PropertyTable"]["Properties"]
    except KeyError:
        raise


def delay_iterations(
    iterable: Iterable[T], waiting_time: float = 60.0, maxsize: int = 400
) -> Iterator[T]:
    """
    Ограничения на запросы, совершаемые в службу PubChem PUG REST:
    * Не больше 5 запросов в секунду.
    * Не больше 400 запросов в минуту.
    * Суммарное время обработки запросов, отправленных в течение минуты, не должно превышать 300
    секунд.
    :param iterable: Последовательности идентификаторов, полученных из функции ``chunked``.
    :param waiting_time: 60 секунд.
    :param maxsize: 400 запросов.
    :return: Генерирует последовательность в соответствии с ограничениями Pubchem.
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


def execute_requests(start: int, stop: int, maxsize: int = 100) -> Iterator[Dict[str, Any]]:
    """
    Извлекает и сохраняет информацию из базы данных Pubchem – 'Compound'; генерирует списки
    идентификаторов равной длины, исходя из параметров в сигнатуре функции, и передает их вместе с
    остальными параметрами, определенными внутри функции, в ``url_builder``, после чего посыпает
    запрос по сгенерированному URL-адресу, учитывая ограничения на количество запросов к серверам
    Pubchem.
    :param start: Первое значение из запрашиваемых CID.
    :param stop: Последнее значение из запрашиваемых CID.
    :param maxsize: Максимальное число идентификаторов в одном запросе, по умолчанию равно 100.
    :return: Генератор запросов к базе данных Pubchem - 'Compound'.
    """
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
            logger.debug("Делаю запрос по URL: {}", url)
            res = request_data_json(url)
        except requests.exceptions.HTTPError:
            logger.error("Ошибка при выполнении запроса.", exc_info=True)
            break
        except KeyError:
            logger.error("Неверный формат ответа от сервера.")
            break
        except BadDomainError:
            logger.error("В данной версии сервиса поиск возможен только по базе данных 'Compound'")
        except BadNamespaceError:
            logger.error("Пространство имен поиска по базе данных 'Compound' задано некорректно")
        except BadOperationError:
            logger.error("Ошибка при составлении операции")
        else:
            logger.debug("Пришел ответ: {}", res)
            for rec in res:
                yield rec
