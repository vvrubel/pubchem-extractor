import pytest

from molecad.builder import *


EXAMPLE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"


def test_input_spec():
    domain = "compound"
    namespace = "cid"
    identifiers = 2244
    res = input_specification(domain, namespace, identifiers)
    assert res == "compound/cid/2244"


def test_build_url():
    input_spec = "compound/cid/2244"
    operation_spec = "property/MolecularFormula,InChIKey"
    res = build_url(input_spec, operation_spec)
    assert res == EXAMPLE


def test_request_data():
    url = EXAMPLE
    params = {}
    res = request_data(url, **params)
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


@pytest.mark.parametrize("inp, expect", [
    ([1, 2, 3], "1,2,3"),
    ([1], "1"),
    ("1", "1"),
])
def test_build_identifiers(inp, expect):
    res = join_w_comma(*inp)
    assert res == expect
