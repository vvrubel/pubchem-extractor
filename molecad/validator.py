from typing import Optional, Sequence

from molecad.types_ import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    PropertyTags,
    SearchPrefix,
    SearchSuffix,
)


def is_not_compound(domain: str) -> bool:
    """
    В текущей версии сервиса доступен запрос свойств молекул из базы данных ``Compound``.
    :param domain: принимает значение ``Domain.COMPOUND``
    :return: если поступивший аргумент равен указанному атрибуту класса, то возвращается ``True``,
    иначе ``False``.
    """
    return domain != Domain.COMPOUND


def is_simple_namespace(prefix: str, suffix: Optional[str]) -> bool:
    """
    Проверяет, что пространство имен поиска состоит только из префикса, а его значение
    удовлетваряет возможным из класса ``NamespCmpd``.
    :param prefix: значение должно принадлежать классу ``NamespCmpd``.
    :param suffix: должен быть равен None.
    :return: если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return suffix is None and isinstance(prefix, NamespCmpd)


def is_namespace_search(prefix: str, suffix: Optional[str]) -> bool:
    """

    :param prefix: значение должно принадлежать классу ``PrefixSearch``.
    :param suffix: значение должно принадлежать классу ``SuffixSearch``.
    :return: если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(suffix, SearchSuffix) and isinstance(prefix, SearchPrefix)


def is_simple_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    """

    :param operation: значение должно принадлежать классу ``Operation``.
    :param tags: должен быть равен None.
    :return: если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(operation, Operation) and tags is None


def is_complex_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    """

    :param operation: значение должно принадлежать классу ``OperationComplex``.
    :param tags: значение должно принадлежать классу ``PropertyTags``.
    :return: если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    return isinstance(operation, OperationComplex) and check_tags(tags)


def check_tags(tags) -> bool:
    """
    Проверяет значение каждого тега на принадлежность к классу ``PropertyTags``.
    :param tags: последовательность, состоящая из строковых значений.
    :return: если условие выполнено, то возвращается ``True``, иначе ``False``.
    """
    checker = []
    if tags is not None:
        checker.append(True)
        for tag in tags:
            checker.append(isinstance(tag, PropertyTags))
    else:
        return False
    return all(checker)
