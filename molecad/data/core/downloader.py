import time
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, TypeVar, Union

import requests
from loguru import logger

from molecad.data.core.utils import (
    chunked,
    generate_ids,
    join_w_comma,
)
from molecad.errors import (
    BadDomainError,
    BadNamespaceError,
    BadOperationError,
)
from molecad.types_ import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    Out,
    PropertyTags,
    SearchPrefix,
    SearchSuffix,
)
from molecad.validator import (
    is_complex_operation,
    is_namespace_search,
    is_not_compound,
    is_simple_namespace,
    is_simple_operation,
)

IdT = TypeVar("IdT", int, str)
T = TypeVar("T")
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def input_specification(domain: str, namespace: str, identifiers: str) -> str:
    """
    Функция, которая форматирует первую часть URL-адреса для запроса в базы данных Pubchem. Эта
    часть определяет, какие записи следует использовать в качестве темы для запроса. Функция
    форматирует параметры для дальнейшей передачи в функцию ``build_url``.
    :param domain: принимает значения из класса ``Domain``.
    .. note:: В текущей версии сервиса доступен поиск только по базе данных ``Compound``.
    :param namespace: принимает значения в зависимости от определенного ранее параметра
    ``domain``. В общем случае может иметь простой и составной характер, поэтому значение должно
    быть предварительно обработано функцией  ``joined_namespaces``.
    :param identifiers: передаваемое значение должно быть целым числом или строкой; если необходимо
    указать последовательность идентификаторов, то они должны быть предварительно переданы в
    функцию ``joined_identifiers``.
    :return: строка, отформатированная по типу "<domain>/<namespace>/<identifiers>",
    которая является первой частью URL-адреса.
    """
    return f"{domain}/{namespace}/{identifiers}"


def operation_specification(operation: str, tags: Optional[str] = None) -> str:
    """
    Часть URL-адреса указывает службе Pubchem, что делать с данными, определенными в первой
    части URL-адреса, например, вы можете получить конкретные свойства соединений. Функция
    форматирует параметры для дальнейшей передачи в функцию ``build_url``.
    :param operation:  принимает значения в зависимости от определенного ранее параметра
    ``domain``. Если значение не определено, то по умолчанию извлекается вся запись.
    :param tags: определяется только в случае если ``operation`` == ``Operation.PROPERTY``,
    иначе не указывается; строка содержит перечень из запрашиваемых тегов, доступных для данной
    операции, которые предварительно были переданы в функцию ``join_w_comma``.
    :return: строка, которая является второй частью URL-адреса.
    """
    if tags is None:
        return f"{operation}"
    else:
        return f"{operation}/{tags}"


def build_url(input_spec: str, operation_spec: str, output_format: str) -> str:
    """
    Команда генерирует строку URL-адреса, необходимую для выполнения запроса в службу PubChem.
    Структура URL-адреса состоит из трех частей, которые предварительно были приведены к
    надлежащему формату.
    :param input_spec: значение полученное из функции ``input_specification``.
    :param operation_spec: значение полученное из функции ``operation_specification``.
    :param output_format: формат выходных данных, принимающий значение из класса ``Out``.
    :return: URL-адрес, используемый в качестве аргумента при вызове функции ``request_data``,
    которая осуществляет запрос в базу данных Pubchem.
    """
    return f"{BASE_URL}/{input_spec}/{operation_spec}/{output_format}"


def joined_namespace(prefix: str, suffix: Optional[str] = None) -> str:
    """
    Функция форматирует входной параметр ``namespace`` для передачи его в функцию
    ``input_specification``. В текущей версии сервиса реализована только для значений из базы
    данных ``Compound``.
    :param prefix: если переменная имеет составной характер, то значение должно быть определено
    из класса ``PrefixSearch``, иначе - из класса ``NamespCmpd``.
    :param suffix: если переменная имеет составной характер, то значение должно быть определено
    из класса ``SuffixSearch``, иначе оно не указывается.
    :return: отформатированное значение ``namespace``, которое в виде строки передается в вызов
    функции ``input_specification``, а при вводе неподходящих параметров кидает ошибку
    ``NamespaceError``.
    """
    if is_simple_namespace(prefix, suffix):
        return f"{prefix}"
    elif is_namespace_search(prefix, suffix):
        return f"{prefix}/{suffix}"
    else:
        raise BadNamespaceError


def prepare_request(
    identifiers: Sequence[IdT],
    domain: Domain = Domain.COMPOUND,
    namespace_prefix: Union[NamespCmpd, SearchPrefix] = NamespCmpd.CID,
    namespace_suffix: Optional[SearchSuffix] = None,
    operation: Union[Operation, OperationComplex] = Operation.RECORD,
    tags: Optional[Sequence[str]] = None,
    output: Out = Out.JSON,
) -> str:
    """
    Подготавливает аргументы и передает их в функцию ``build_url``.
    .. note:: В текущей версии сервиса доступен запрос свойств молекул из базы данных ``Compound``.
    :param identifiers: число, строка или последовательность из строковых значений/чисел,
    которые интерпретируются как идентификаторы соединений; последовательность идентификаторов может
    быть получена в результате выполнения функции ``generate_ids`` или ``chunked``, полученные
    значения передаются в функцию ``join_w_comma``, а затем в ``input_specification``.
    :param domain: принимает значение из класса ``Domain`` и напрямую передается в функцию
    ``input_specification``. По умолчанию установлено значение ``Domain.COMPOUND``.
    .. note:: В текущей версии сервиса доступен поиск только по базе данных ``Compound``.
    :param namespace_prefix: передается в функцию ``joined_namespace`` в качестве ``prefix``,
    который может принимать значения из классов ``NamespCmpd`` и ``PrefixSearch`` в данной
    версии сервиса. Значение по умолчанию - ``NamespCmpd.CID``.
    .. note:: В текущей версии сервиса доступен поиск только по пространству имен ``CID``.
    :param namespace_suffix: передается в функцию ``joined_namespace`` в качестве ``suffix`` и
    принимает значения из класса ``SuffixSearch``; по умолчанию - None
    .. note:: В текущей версии сервиса значение параметра - None.
    :param operation: принимает значения в зависимости от определенного ранее параметра
    ``domain``. Аргумент может принимать значения из классов ``Operation``, ``OperationComplex``.
    Если значение не определено, то по умолчанию извлекается вся запись, т.е. значение равно
    ``Operation.RECORD``.
    .. note:: В текущей версии сервиса доступна операции, поиск по которым выполняется в базе
    данных ``Compound``, а формат их выходных данных должен представлять собой JSON.
    :param tags: определяется только в случае если ``operation`` == ``Operation.PROPERTY``,
    иначе не указывается; строка содержит перечень из запрашиваемых тегов, доступных для данной
    операции, которые передаются в функцию ``join_w_comma`` и далее в ``operation_specification``.
    :param output: принимает значение из класса ``Out`` и напрямую передается в функцию
    ``build_url``; по умолчанию - ``Out.JSON``.
    .. note:: В текущей версии сервиса получение ответа от сервера возможно только в формате JSON.
    :return: URL-адрес, который используется для запроса в базу данных Pubchem.
    """
    if is_not_compound(domain):
        raise BadDomainError
    joined_identifiers = join_w_comma(*identifiers)
    namespace = joined_namespace(namespace_prefix, namespace_suffix)
    input_spec = input_specification(domain, namespace, joined_identifiers)

    if is_complex_operation(operation, tags) and tags is not None:
        joined_tags = join_w_comma(*tags)
        operation_spec = operation_specification(operation, joined_tags)
    elif is_simple_operation(operation, tags):
        operation_spec = operation_specification(operation)
    else:
        raise BadOperationError
    url = build_url(input_spec, operation_spec, output)
    return url


def request_data_json(url: str, **params: str) -> List[Dict[str, Any]]:
    """
    Функция отправляет синхронный запрос "GET" к серверам PUG REST баз данных Pubchem.
    :param url: возвращается из функции ``prepare_request`` и является обязательным для всех
    запросов.
    :param params: может быть пустым словарем, если параметры ``operation`` не требуют иного,
    иначе передается словарь со строковыми значениями, которые интерпретируются как ``operation
    options``. Если бы вы создавали URL-адрес вручную, пары значений из словаря
    записыварись бы в конец URL-адреса после знака ``?`` в виде ``key=value``, а при наличии
    более одной пары объединялись знаком ``&``.
    :return: ответ от сервера .. note:: В текущей версии сервиса реализуется получения ответа
    только в формате JSON .
    """
    response = requests.get(url, params=params).json()
    return response["PropertyTable"]["Properties"]


def delay_iterations(
    iterable: Iterable[T], waiting_time: float = 60.0, maxsize: int = 400
) -> Iterator[T]:
    """
    Ограничения на запросы, совершаемые в службу PubChem REST:
    * Не больше 5 запросов в секунду.
    * Не больше 400 запросов в минуту.
    * Суммарное время обработки запросов, отправленных в течение минуты, не должно превышать 300
    секунд.
    :param iterable: последовательности идентификаторов, полученных из функции ``generate_ids``
    или ``chunked``.
    :param waiting_time: 60 секунд.
    :param maxsize: 400 запросов.
    :return: генерирует последовательность в соответствии с ограничениями серверов Pubchem.
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


def execute_request(start: int, stop: int, maxsize: int) -> Dict[Any, Dict[str, Any]]:
    """
    .. note:: В текущей версии сервиса доступен запрос свойств молекул из базы данных ``Compound``.
    Аргументы функции ``generate_ids(start, stop)`` по умолчанию равны 1 и 201 соответственно и
    могут не указываться явно, что соответствует тестовым запросам к серверу Pubchem;
    в случае формирования базы данных эти значения должны быть явно указаны в качестве
    аргументов, как 1 и 500001 соответственно. Для последующих запросов ``start`` = ``stop`` от
    предыдущего запроса, а ``stop`` увеличивается на 500000.
    Второй аргумент, передаваемый в функцию ``chunked`` - ``chunk_size``, в рамках для тестовых
    запросов по умолчанию равен 100 и может не указываться явно; при формировании базы данных со
    свойствами молекул должен быть равен 1000.
    """
    domain = Domain.COMPOUND
    namespace_prefix = NamespCmpd.CID
    namespace_suffix = None
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

    data = {}
    t_start = time.monotonic()
    chunks = chunked(generate_ids(start, stop), maxsize)
    for i in delay_iterations(chunks):
        url = prepare_request(
            i, domain, namespace_prefix, namespace_suffix, operation, tags, output
        )
        logger.debug("Делаю запрос по URL: {}", url)
        try:
            res = request_data_json(url)
        except requests.HTTPError:
            logger.error("Ошибочка вышла: {}", exc_info=True)
            break
        else:
            logger.debug("Пришел ответ: {}", res)
            for k, v in zip(i, res):
                data[k] = v
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    logger.info("Время, затраченное на операцию, равно {}", t_run)

    return data
