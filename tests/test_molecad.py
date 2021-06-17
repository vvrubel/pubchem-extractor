import pytest

from molecad.builder import (
    build_url,
    input_specification,
    operation_specification
)
from molecad.sample_sync_requests import (
    chunked,
    generate_ids,
    request_data_json
)
from molecad.main import (
    join_w_comma,
    prepare_request
)
from molecad import __version__

EXAMPLE1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"
EXAMPLE2 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG"
EXAMPLE3 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
EXAMPLE4 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
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


def test_input_spec():
    domain = "compound"
    namespace = "cid"
    identifiers = 2244
    res = input_specification(domain, namespace, identifiers)
    assert res == "compound/cid/2244"


def test_input_spec_ids():
    domain = "compound"
    namespace = "cid"
    identifiers = "1,2,3"
    res = input_specification(domain, namespace, identifiers)
    assert res == "compound/cid/1,2,3"


def test_operation_spec():
    operation = "record"
    res = operation_specification(operation)
    assert res == "record"


def test_operation_spec_w_tags():
    operation = "property"
    tags = "MolecularFormula,InChIKey"
    res = operation_specification(operation, tags)
    assert res == "property/MolecularFormula,InChIKey"


def test_build_url():
    input_spec = "compound/cid/2244"
    operation_spec = "record"
    output = "PNG"
    res = build_url(input_spec, operation_spec, output)
    assert res == EXAMPLE2


def test_build_url_w_tags():
    input_spec = "compound/cid/2244"
    operation_spec = "property/MolecularFormula,InChIKey"
    output = "JSON"
    res = build_url(input_spec, operation_spec, output)
    assert res == EXAMPLE1


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
    ids.extend(generate_ids())
    expectation = [int(i) for i in range(1, 11)]
    assert ids == expectation


def test_chunked():
    ids = generate_ids()
    chunks = []
    chunks.extend(chunked(ids, 2))
    assert chunks == [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]


def test_prepare_request():
    domain = "compound"
    namespace = "cid"
    identifiers = (1,)
    operation = "property"
    output = "JSON"
    url = prepare_request(domain, namespace, identifiers, operation, output)
    assert url == BAD_EXAMPLE


def test_prepare_request_w_tags():
    domain = "compound"
    namespace = "cid"
    identifiers = (1,)
    operation = "property"
    tags = ["MolecularFormula", "MolecularWeight", "IUPACName", "CanonicalSMILES"]
    output = "JSON"
    url = prepare_request(domain, namespace, identifiers, operation, output, tags)
    assert url == EXAMPLE3


def test_prepare_request_chunk():
    domain = "compound"
    namespace = "cid"
    identifiers = []
    identifiers.extend(chunked(generate_ids(), 2))
    operation = "property"
    tags = ["MolecularFormula", "MolecularWeight", "IUPACName", "CanonicalSMILES"]
    output = "JSON"
    url = prepare_request(domain, namespace, identifiers[0], operation, output, tags)
    assert url == EXAMPLE4
