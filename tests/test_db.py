import pytest

from molecad.settings import setup
from molecad.utils import index_name_extractor
from molecad.validator import (
    is_index_exist,
)


def test_is_index_exist():
    collection, schema = setup.get_prod_collection
    res = is_index_exist(collection, "CID")
    assert res is True


def test_is_index_not_exist():
    test = setup.get_dev_collection
    res = is_index_exist(test, "CID")
    assert res is False


@pytest.mark.parametrize(
    "inp, expect",
    [
        ([("CID", 1)], ["CID_1"]),
        ([("CID", -1)], ["CID_-1"]),
    ],
)
def test_index_name_extractor(inp, expect):
    res = index_name_extractor(inp)
    assert res == expect
