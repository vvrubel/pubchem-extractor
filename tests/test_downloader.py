import pytest

from molecad.downloader import chunked, generate_ids, request_data_json
from molecad.types_ import Domain, NamespCmpd, Operation, OperationComplex, PropertyTags
from molecad.utils import concat
from molecad.validator import check_tags, is_complex_operation, is_simple_operation

EXAMPLE1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"


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
        ((1, "two", 3), "1,two,3"),
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
        ("compound", "cid", 2244, "compound/cid/2244"),
        ("compound", "cid", "1,2,3", "compound/cid/1,2,3"),
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


@pytest.mark.parametrize(
    "tags, expect",
    [
        (None, False),
        ((PropertyTags.MOLECULAR_FORMULA, PropertyTags.CANONICAL_SMILES), True),
    ],
)
def test_check_tags(tags, expect):
    res = check_tags(tags)
    assert res == expect


def test_request_data_json():
    url = EXAMPLE1
    params = {}
    res = request_data_json(url, **params)
    expectation = [
        {
            "CID": 2244,
            "MolecularFormula": "C9H8O4",
            "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        }
    ]
    assert res == expectation

