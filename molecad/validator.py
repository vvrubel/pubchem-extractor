from molecad.data.request import *
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
    return domain != Domain.COMPOUND


def is_simple_namespace(prefix: str, suffix: Optional[str]) -> bool:
    return suffix is None and isinstance(prefix, NamespCmpd)


def is_namespace_search(prefix: str, suffix: Optional[str]) -> bool:
    return isinstance(suffix, SearchSuffix) and isinstance(prefix, SearchPrefix)


def is_simple_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    return isinstance(operation, Operation) and tags is None


def is_complex_operation(operation: str, tags: Optional[Sequence[str]]) -> bool:
    return isinstance(operation, OperationComplex) and tags is not None


def check_tags(tags: Sequence[str]) -> bool:
    """
    Проверяет значение каждого тега на принадлежность к классу ``PropertyTags``.
    :param tags: последовательность строковых значений.
    :return: если все элементы принадлежат к указанному классу, то возвращает ``True``,
    иначе ``False``.
    """
    checker = []
    for tag in tags:
        checker.append(isinstance(tag, PropertyTags))
    return all(checker)





# def has_allowed_extension(filename: str, allowed_extensions=("smi",)) -> bool:
#     if filename and all(filename.endswith(allow_ext) for allow_ext in allowed_extensions):
#         return True
#     raise WrongExtensionError
#
#
# def has_content(content: bytes) -> bool:
#     if content:
#         return True
#     raise EmptyFileError
#
#
# def without_empty_strings(line: bytes, line_number: int) -> bool:
#     if not line:
#         raise ContainsEmptyStringsError(line=line_number)
#     return True
#
#
# def has_correct_columns_number(line: bytes, allowed_columns_number: int = 2) -> bool:
#     columns = line.split(b" ")
#     if len(columns) != allowed_columns_number:
#         raise IncorrectNumberOfColumnsError(allowed_columns_number=allowed_columns_number)
#     return True
#
#
# def has_valid_molecule(line: bytes) -> bool:
#     smiles = line.split(b" ")[0]
#     molecule = Chem.MolFromSmiles(smiles)
#     if molecule is None:
#         raise InvalidMoleculeError(wrong_molecule=smiles.decode())
#     return True
#
#
# def each_line_is_valid(content: bytes) -> bool:
#     lines: list = content.strip().splitlines()
#
#     for line_number, line in enumerate(lines):
#         without_empty_strings(line, line_number + 1)
#         has_correct_columns_number(line)
#         has_valid_molecule(line)
#
#     return True
