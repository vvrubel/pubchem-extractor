import pytest

from molecad.data.downloader import (
    chunked,
    generate_ids,
    request_property_data_json,
)
from molecad.data.utils import (
    concat,
)
from molecad.types_ import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    PropertyTags,
)
from molecad.validator import (
    is_complex_operation,
    is_simple_operation,
)

EXAMPLE1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"
EXAMPLE2 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG"
EXAMPLE3 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE4 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"


@pytest.mark.parametrize(
    "start, stop, expect",
    [
        (1, 11, [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]),
        (1, 10, [[1, 2], [3, 4], [5, 6], [7, 8], [9]]),
    ],
)
def test_chunked(start, stop, expect):
    ids = generate_ids(start, stop)
    chunks = []
    chunks.extend(chunked(ids, 2))
    assert chunks == expect


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1,2,3"),
        ((1, 2, 3), "1,2,3"),
        ((1, 2, 3,), "1,2,3"),
        ([1], "1"),
        ("1", "1"),
        ([], ""),
    ],
)
def test_concat_w_comma(inp, expect):
    res = concat(*inp, sep=",")
    assert res == expect


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1/2/3"),
        ([1], "1"),
    ],
)
def test_concat_w_slash(inp, expect):
    res = concat(*inp)
    assert res == expect


@pytest.mark.parametrize(
    "domain, namespace, ids, expect",
    [
        (Domain.COMPOUND, NamespCmpd.CID, 2244, "compound/cid/2244"),
        (Domain.COMPOUND, NamespCmpd.CID, "1,2,3", "compound/cid/1,2,3"),
    ],
)
def test_input_specification(domain, namespace, ids, expect):
    res = concat(domain, namespace, ids)
    assert res == expect


def test_simple_operation():
    operation = Operation.RECORD
    tags = None
    expect = "record"
    if is_simple_operation(operation, tags):
        res = operation
        assert res == expect


def test_complex_operation():
    operation = OperationComplex.PROPERTY
    tags = (PropertyTags.MOLECULAR_FORMULA, PropertyTags.INCHI)
    expect = "property/MolecularFormula,InChI"
    if is_complex_operation(operation, tags):
        joined_tags = concat(*tags, sep=",")
        res = concat(operation, joined_tags)
        assert res == expect


def test_request_data_json():
    url = EXAMPLE1
    params = {}
    res = request_property_data_json(url, **params)
    expectation = [
        {
            "CID": 2244,
            "MolecularFormula": "C9H8O4",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        }
    ]
    assert res == expectation



#
#
# def test_prepare_request_w_tags():
#     domain = Domain.COMPOUND
#     namespace_prefix = NamespCmpd.CID
#     namespace_suffix = None
#     i = "1"
#     operation = OperationComplex.PROPERTY
#     tags = (
#         PropertyTags.MOLECULAR_FORMULA,
#         PropertyTags.MOLECULAR_WEIGHT,
#         PropertyTags.IUPAC_NAME,
#         PropertyTags.CANONICAL_SMILES,
#     )
#     output = Out.JSON
#     url = url_builder(i, domain, namespace_prefix, namespace_suffix, operation, tags, output)
#     assert url == EXAMPLE3

#
# def test_prepare_request_chunk():
#     domain = Domain.COMPOUND
#     namespace_prefix = NamespCmpd.CID
#     namespace_suffix = None
#     ids = generate_ids(1, 3)
#     operation = OperationComplex.PROPERTY
#     tags = (
#         PropertyTags.MOLECULAR_FORMULA,
#         PropertyTags.MOLECULAR_WEIGHT,
#         PropertyTags.IUPAC_NAME,
#         PropertyTags.CANONICAL_SMILES,
#     )
#     output = Out.JSON
#     for chunk in chunked(ids, 2):
#         identifiers = chunk
#         url = url_builder(
#             identifiers, domain, namespace_prefix, namespace_suffix, operation, tags, output
#         )
#         assert url == EXAMPLE4
