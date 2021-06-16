import argparse
import sys
from molecad.sync_requests import request_data
from molecad.builder import *


EXAMPLE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"


def build_parser():
    parser = argparse.ArgumentParser(description="This module can help you to specify UPL parameters")
    parser.add_argument('--async', action=argparse.BooleanOptionalAction)
    # parser.add_argument('-d', '--domain', choices=["compound", "substance", "assay"], default="compound")
    # parser.add_argument('-n', '--namespace', choices=["cid", "name", "smiles"], default="cid")
    parser.add_argument('--ids', type=int, nargs='*')
    return parser


def build_request(identifiers, tags):
    joined_identifiers = join_w_comma(*identifiers)
    input_spec = input_specification(joined_identifiers)
    joined_tags = join_w_comma(*tags)
    operation_spec = operation_specification(joined_tags)
    url = build_url(input_spec, operation_spec,)
    return url


def execute_request(url, params):
    print(url)
    res = request_data(url, **params)
    print(res)
    return res


def main(*args):
    parser = build_parser()
    params = parser.parse_args(args)
    print(params)
    # ids = params.ids # generate_ids()
    # chunks = chunked(ids, 2)
    tags = ["MolecularFormula", "MolecularWeight", "IUPACName", "CanonicalSMILES"]
    # additional_tags = ["InChI", "XLogP", "TPSA", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "Volume3D"
    results = {}
    for id in ids:
        url = build_request([id], tags)
        res = execute_request(url, {})
        results[id] = res


if __name__ == '__main__':
    main(*sys.argv[1:])
