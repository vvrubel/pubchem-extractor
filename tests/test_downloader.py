import pytest

from molecad.data.downloader import (
    build_url,
    chunked,
    generate_ids,
    input_specification,
    operation_specification,
    prepare_request,
    request_data_json,
)
from molecad.data.utils import (
    concat,
    join_w_comma,
)
from molecad.types_ import (
    Domain,
    NamespCmpd,
    Operation,
    OperationComplex,
    Out,
    PropertyTags,
)

EXAMPLE1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"
EXAMPLE2 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG"
EXAMPLE3 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE4 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE5 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG?image_size=large"
BAD_EXAMPLE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/JSON"


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1,2,3"),
        ([1], "1"),
        ("1", "1"),
    ],
)
def test_join_w_comma(inp, expect):
    res = join_w_comma(*inp)
    assert res == expect


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1,2,3"),
        ([1], "1"),
        ("1", "1"),
    ],
)
def test_concat_w_comma(inp, expect):
    res = concat(*inp, sep=",")
    assert res == expect


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1–2–3"),
        ([1, 2], "1–2"),
        ([1], "1"),
        ("1", "1"),
    ],
)
def test_concat_w_dash(inp, expect):
    res = concat(*inp, sep="–")
    assert res == expect


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([1, 2, 3], "1/2/3"),
        ([1], "1"),
        ("1", "1"),
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
        ("compound", NamespCmpd.CID, "1,2,3", "compound/cid/1,2,3"),
        ("compound", "cid", "1,2,3", "compound/cid/1,2,3"),
        (
            Domain.COMPOUND,
            "xref/PatentID",
            "EP0001699A2",
            "compound/xref/PatentID/EP0001699A2",
        ),
    ],
)
def test_input_specification(domain, namespace, ids, expect):
    res = input_specification(domain, namespace, ids)
    assert res == expect


@pytest.mark.parametrize(
    "op, tags, expect",
    [
        ("record", None, "record"),
        ("property", "MolecularFormula,InChIKey", "property/MolecularFormula,InChIKey"),
    ],
)
def test_operation_specification(op, tags, expect):
    res = operation_specification(op, tags)
    assert res == expect


@pytest.mark.parametrize(
    "inp, op, out, expect",
    [
        ("compound/cid/2244", "property/MolecularFormula,InChIKey", "JSON", EXAMPLE1),
        ("compound/cid/2244", "record", "PNG", EXAMPLE2),
        ("compound/cid/2244", Operation.RECORD, Out.PNG, EXAMPLE2),
    ],
)
def test_build_url(inp, op, out, expect):
    res = build_url(inp, op, out)
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


def test_generate_ids():
    ids = []
    ids.extend(generate_ids(1, 11))
    expectation = [int(i) for i in range(1, 11)]
    assert ids == expectation


def test_chunked():
    ids = generate_ids(1, 11)
    chunks = []
    chunks.extend(chunked(ids, 2))
    assert chunks == [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]


def test_prepare_request_w_tags():
    domain = Domain.COMPOUND
    namespace_prefix = NamespCmpd.CID
    namespace_suffix = None
    i = "1"
    operation = OperationComplex.PROPERTY
    tags = (
        PropertyTags.MOLECULAR_FORMULA,
        PropertyTags.MOLECULAR_WEIGHT,
        PropertyTags.IUPAC_NAME,
        PropertyTags.CANONICAL_SMILES,
    )
    output = Out.JSON
    url = prepare_request(i, domain, namespace_prefix, namespace_suffix, operation, tags, output)
    assert url == EXAMPLE3


def test_prepare_request_chunk():
    domain = Domain.COMPOUND
    namespace_prefix = NamespCmpd.CID
    namespace_suffix = None
    ids = generate_ids(1, 3)
    operation = OperationComplex.PROPERTY
    tags = (
        PropertyTags.MOLECULAR_FORMULA,
        PropertyTags.MOLECULAR_WEIGHT,
        PropertyTags.IUPAC_NAME,
        PropertyTags.CANONICAL_SMILES,
    )
    output = Out.JSON
    for chunk in chunked(ids, 2):
        identifiers = chunk
        url = prepare_request(
            identifiers, domain, namespace_prefix, namespace_suffix, operation, tags, output
        )
        assert url == EXAMPLE4
