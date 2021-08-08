from typing import Optional, Sequence

from molecad.cli_tools.url_parts import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    PropertyTags,
    SearchPrefix,
    SearchSuffix,
)


def is_compound(domain: str) -> bool:
    return domain == Domain.COMPOUND


def is_simple_namespace(prefix: str, suffix: Optional[str]) -> bool:
    return suffix is None and isinstance(prefix, NamespCmpd)


def is_namespace_search(prefix: str, suffix: Optional[str]) -> bool:
    return isinstance(suffix, SearchSuffix) and isinstance(prefix, SearchPrefix)


def is_simple_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    return isinstance(operation, Operation) and tags is None


def is_complex_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    return isinstance(operation, OperationComplex) and check_tags(tags)


def check_tags(tags: Optional[Sequence[str]]) -> bool:
    if tags is None:
        return False
    else:
        return all(isinstance(tag, PropertyTags) for tag in tags)


# TODO: add validator on query params
