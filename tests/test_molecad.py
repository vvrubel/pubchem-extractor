import pytest

from molecad.sample_sync_requests import *
from molecad import __version__


def test_version():
    assert __version__ == '0.1.0'


class TestRequestData:
    def test_request_data(self):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/JSON"
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

    # def test_request_data_w_params(self):
    #     url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/record/PNG"
    #     params = {'image_size': 'small'}
    #     res = request_data(url, **params)
    #     assert res.headers["Content-Type"] == "image/png"
    #     assert res.status_code == 200
