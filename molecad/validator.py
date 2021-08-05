from typing import Optional, Sequence

from molecad.cli_tools.downloader_types import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    PropertyTags,
    SearchPrefix,
    SearchSuffix,
)


def is_compound(domain: str) -> bool:
    """
    В текущей версии сервиса доступен запрос свойств молекул из базы данных "Compound".
    :param domain: Значение должно быть равным ``Domain.COMPOUND``
    :return: Если поступивший аргумент равен ``Domain.COMPOUND``, то возвращается ``True``,
    иначе ``False``.
    """
    return domain == Domain.COMPOUND


def is_simple_namespace(prefix: str, suffix: Optional[str]) -> bool:
    """
    Проверяет, что пространство имен поиска состоит только из префикса, а его значение принадлежит
    классу ``NamespCmpd``.
    :param prefix: Значение должно принадлежать классу ``NamespCmpd``.
    :param suffix: Должен быть равен ``None``.
    :return: Если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return suffix is None and isinstance(prefix, NamespCmpd)


def is_namespace_search(prefix: str, suffix: Optional[str]) -> bool:
    """
    Проверяет, что пространство имен поиска составлено корректно.
    :param prefix: Значение должно принадлежать классу ``PrefixSearch``.
    :param suffix: Значение должно принадлежать классу ``SuffixSearch``.
    :return: Если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(suffix, SearchSuffix) and isinstance(prefix, SearchPrefix)


def is_simple_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    """
    Проверяет,что операция является простой и составлена корректно.
    :param operation: Значение должно принадлежать классу ``Operation``.
    :param tags: Должен быть равен ``None``.
    :return: Если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(operation, Operation) and tags is None


def is_complex_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    """
    Проверяет,что операция является составной и сочетание аргументов корректно.
    :param operation: Значение должно принадлежать классу ``OperationComplex``.
    :param tags: Итерируемый объект отправляется на проверку в функцию ``check_tags``.
    :return: Если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(operation, OperationComplex) and check_tags(tags)


def check_tags(tags: Optional[Sequence[str]]) -> bool:
    """
    Проверяет значение каждого тега на принадлежность к классу ``PropertyTags``.
    :param tags: Последовательность, состоящая из строковых значений.
    :return: Если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    if tags is None:
        return False
    else:
        return all(isinstance(tag, PropertyTags) for tag in tags)
