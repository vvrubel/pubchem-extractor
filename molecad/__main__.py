import argparse
import requests
import sys
import time
from typing import (
    Any,
    Iterable,
    Optional,
    TypeVar
)

from molecad.builder import (
    build_url,
    input_specification,
    operation_specification
)
from molecad.sample_sync_requests import (
    generate_ids,
    chunked,
    delay_iterations,
    request_data_json
)

T = TypeVar("T")
EXAMPLE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"


def join_w_comma(*args: object) -> str:
    """
    Formats all iterable arguments to a string with comma.
    :param args: object to indicate that a value could be any type in a typesafe manner
    :return: string of comma-separated values without whitespaces
    """
    return ",".join(map(str, args))


def prepare_request(
        domain: str,
        namespace: str,
        identifiers: Iterable[T],
        operation: str,
        output: str,
        property_tags: Optional[list] = None
) -> str:
    """
    Prepares all arguments to be passed to URL builder. Mainly, joins all list arguments with comma to string.
    :param domain:
    :param namespace:
    :param identifiers: yielded value from
    :param operation:
    :param output:
    :param property_tags:
    :return:
    """
    joined_identifiers = join_w_comma(*identifiers)
    input_spec = input_specification(domain, namespace, joined_identifiers)
    if property_tags is None:
        operation_spec = operation_specification(operation)
    else:
        joined_tags = join_w_comma(*property_tags)
        operation_spec = operation_specification(operation, joined_tags)
    url = build_url(input_spec, operation_spec, output)
    return url


def execute_request(url, params):
    print(url)
    res = request_data_json(url, **params)
    return res


def main(*args: Any):
    # # parser = build_parser()
    # # params = parser.parse_args(args)
    # # print(params)
    # # ids = params.ids # generate_ids()
    # # chunks = chunked(ids, 2)
    # # out_box = ["JSON", "XML", "SDF", "CSV", "PNG", "TXT"]

    domain = "compound"
    namespace = "cid"
    operation = "property"
    output = "JSON"
    tags = ["MolecularFormula", "MolecularWeight", "IUPACName", "CanonicalSMILES"]
    # additional_tags = ["InChI", "XLogP", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "Volume3D"]

    results = {}
    t_start = time.monotonic()
    for i in delay_iterations(generate_ids(), 60.0, 400):
        url = prepare_request(domain, namespace, [i], operation, output, tags)
        try:
            res = execute_request(url, {})
        except requests.HTTPError:
            break
        else:
            results[i] = res
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    print(t_run)


if __name__ == '__main__':
    main(*sys.argv[1:])
