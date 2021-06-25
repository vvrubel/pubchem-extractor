import pytest

from molecad import __version__
from molecad.builder import (
    build_url,
    input_specification,
    join_w_comma,
    joined_namespace,
    operation_specification,
    prepare_request,
)
from molecad.downloader import *
from molecad.types_ import (
    IdT,
    Domain,
    NamespaceComp,
    Operations,
    Out,
    PropertyTags,
    SearchSuffix,
    WrongValue,
    Xref
)

EXAMPLE1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"
EXAMPLE2 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG"
EXAMPLE3 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE4 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE5 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG?image_size=large"
BAD_EXAMPLE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/JSON"


def test_version():
    assert __version__ == '0.1.0'


@pytest.mark.parametrize("inp, expect", [
    ([1, 2, 3], "1,2,3"),
    ([1], "1"),
    ("1", "1"),
])
def test_join_w_comma(inp, expect):
    res = join_w_comma(*inp)
    assert res == expect


@pytest.mark.parametrize("pref, suf, expect", [
    (NamespaceComp.CID, None, "cid"),
    (NamespaceComp.IDENTITY, SearchSuffix.CID, "identity/cid"),
    (NamespaceComp.XREF, Xref.PATENT_ID, "xref/PatentID"),
])
def test_only_prefix(pref, suf, expect):
    res = joined_namespace(pref, suf)
    assert res == expect


@pytest.mark.parametrize("domain, namespace, ids, expect", [
    (Domain.COMPOUND.value, NamespaceComp.CID.search, 2244, "compound/cid/2244"),
    (Domain.COMPOUND.value, NamespaceComp.CID.search, "1,2,3", "compound/cid/1,2,3"),
    ("compound", NamespaceComp.CID.search, "1,2,3", "compound/cid/1,2,3"),
    ("compound", "cid", "1,2,3", "compound/cid/1,2,3"),
    (Domain.COMPOUND.value, "xref/PatentID", "EP0001699A2", "compound/xref/PatentID/EP0001699A2"),
])
def test_input_specification(domain, namespace, ids, expect):
    res = input_specification(domain, namespace, ids)
    assert res == expect


@pytest.mark.parametrize("op, tags, expect", [
    ("record", None, "record"),
    ("property", "MolecularFormula,InChIKey", "property/MolecularFormula,InChIKey"),
])
def test_operation_specification(op, tags, expect):
    res = operation_specification(op, tags)
    assert res == expect


@pytest.mark.parametrize("inp, op, out, expect", [
    ("compound/cid/2244", "property/MolecularFormula,InChIKey", "JSON", EXAMPLE1),
    ("compound/cid/2244", "record", "PNG", EXAMPLE2),
    ("compound/cid/2244", Operations.RECORD.action, Out.PNG.value, EXAMPLE2),
])
def test_build_url(inp, op, out, expect):
    res = build_url(inp, op, out)
    assert res == expect


def test_request_data_json():
    url = EXAMPLE1
    params = {}
    res = request_data_json(url, **params)
    expectation = {
        'PropertyTable': {
            'Properties': [
                {
                    'CID': 2244,
                    'MolecularFormula': 'C9H8O4',
                    'InChIKey': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
                }
            ]
        }
    }
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


def test_prepare_request():
    domain = "compound"
    namespace = "cid"
    identifiers = [1]
    operation = "property"
    output = "JSON"
    url = prepare_request(domain, namespace, identifiers, operation, output)
    assert url == BAD_EXAMPLE
#
#
# def test_prepare_request_w_tags():
#     domain = Domain.COMPOUND
#     namespace = NamespaceComp.CID
#     identifiers = (1,)
#     operation = Operations.PROPERTY
#     tags = (
#         PropertyTags.MOLECULAR_FORMULA,
#         PropertyTags.MOLECULAR_WEIGHT,
#         PropertyTags.IUPAC_NAME,
#         PropertyTags.CANONICAL_SMILES,
#     )
#     output = Out.JSON
#     url = prepare_request(
#         domain,
#         namespace,
#         identifiers,
#         operation,
#         output,
#         tags)
#     assert url == EXAMPLE3
#
#
# def test_prepare_request_chunk():
#     domain = Domain.COMPOUND
#     namespace = NamespaceComp.CID
#     ids = generate_ids(1, 3)
#     operation = Operations.PROPERTY
#     tags = (
#         PropertyTags.MOLECULAR_FORMULA,
#         PropertyTags.MOLECULAR_WEIGHT,
#         PropertyTags.IUPAC_NAME,
#         PropertyTags.CANONICAL_SMILES,
#     )
#     output = Out.JSON
#     for chunk in chunked(ids, 2):
#         identifiers = chunk
#         url = prepare_request(
#             domain,
#             namespace,
#             identifiers,
#             operation,
#             output,
#             tags)
#         assert url == EXAMPLE4
