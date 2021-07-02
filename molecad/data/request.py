import time
from typing import (Any, Dict, Iterable, Iterator, List, Optional, Sequence, TypeVar, Union)

import requests

from molecad.errors import (
    BadDomainError,
    BadNamespaceError,
    OperationError,
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
    check_tags,
    is_complex_operation,
    is_namespace_search,
    is_not_compound,
    is_simple_namespace,
    is_simple_operation,
)

IdT = TypeVar("IdT", int, str)
T = TypeVar("T")
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def generate_ids(start: int = 1, stop: int = 201) -> Iterator[int]:
    """
    Простой генератор значений CID.
    .. note:: В рамках формирования базы данных интервал идентификаторов был равен (1, 500001).
    :param start: по умолчанию равно 1, но может быть заменено на любое положительное число -
    для скачивания порциями равно значению ``stop`` в предыдущей порции загрузки.
    :param stop: любое значение до 156 миллионов; в тестовом режиме установлено значение 201 для
    получения быстрого результата.
    :return: генератор целых положительных чисел.
    """
    for n in range(start, stop):
        yield n


def chunked(iterable: Iterable[T], maxsize: int = 100) -> Iterator[List[T]]:
    """
    Принимает итерируемый объект со значениями одинакового типа и делит его на чанки одинкаовой
    длины, равной ``maxsize``.
    :param iterable: последовательность, пришедшая из функции ``generate_ids``.
    :param maxsize: максимальное число элементов в чанке, для формирования базы данных равно
    1000, тестовые запросы должны выполняться со значением 100.
    :return: выбрасывает списоки элементов - чанки, которые затем необходимо передать в функцию
    ``join_w_comma()``, после чего полученная строка может быть передана в ``input_specification``.
    """
    chunk = []
    for i in iterable:
        chunk.append(i)
        if len(chunk) >= maxsize:
            yield chunk
            chunk = []


def delay_iterations(
    iterable: Iterable[T],
    waiting_time: float = 60.0,
    maxsize: int = 400
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
    :return: генерирует последовательность в соотвествии с ограничениями серверов Pubchem.
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


def operation_specification(
    operation: str, tags: str = None
) -> str:
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


def join_w_comma(*args: T) -> str:
    """
    Функция форматирует входящую последовательность аргументов, соединяя элементы запятой без
    пробелов и приводит полученное к строке.
    :param args: любая последовательность аргументов одинакового типа.
    :return: строка разделенных запятой и без пробела значений.
    """
    return ",".join(f"{i}" for i in args)


def joined_namespace(
    prefix: str,
    suffix: Optional[str] = None
) -> str:
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
    tags: Optional[Sequence[PropertyTags]] = None,
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
    ``domain``. Если значение не определено, то по умолчанию извлекается вся запись,
    т.е. значение равно ``Operation.RECORD``.
    .. note:: В текущей версии сервиса доступена операции, поиск по которым выполняется в базе
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
    if is_complex_operation(operation, tags) and check_tags(tags):
        joined_tags = join_w_comma(tags)
        operation_spec = operation_specification(operation, joined_tags)
    elif is_simple_operation(operation, tags):
        operation_spec = operation_specification(operation)
    else:
        raise OperationError
    url = build_url(input_spec, operation_spec, output)
    return url


def request_data_json(url: str, **params: str) -> Dict[str, Any]:
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
    return response
