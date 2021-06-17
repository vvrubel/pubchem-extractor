import argparse
from typing import (
    Optional,
    Union
)


def input_specification(domain: str, namespace: str, identifiers: Union[int, str]) -> str:
    """
    This function builds the first part of URL setting up input identifiers from PubChem databases.
    :param domain: string refers to database - compound | substance | assay
    :param namespace: string can assign values depending on the specified domain.
    #TODO make description
    if domain == `compound`:
    if domain == `substance`:
    if domain == `assay`:
    :param identifiers:
    :return: string is the first part of URL and formatted as `<domain>/<namespace>/<identifiers>`
    """
    return f"{domain}/{namespace}/{identifiers}"


def operation_specification(operation: str, property_tags:  Optional[str] = None) -> str:
    """
    The function dictates what to do with the input records – what information about them you want to retrieve.
    :param operation: string describes action to execute with input identifiers.
    :param property_tags: iterable or string are passed if only `operation == "property"`
    :return: string is the second part of URL
    """
    if property_tags is None:
        return f"{operation}"
    else:
        return f"{operation}/{property_tags}"


def build_url(input_spec: str, operation_spec: str, output_format: str) -> str:
    """
    The command generates the URL string needed to request data from PubChem’s PUG REST service.
    URL structure consists of three pre-build parts assigned as parameters when invoking the function.
    # TODO
    :param input_spec: string represented as "<domain>/<namespace>/<identifiers>"
    :param operation_spec: string
    :param output_format: takes one value from the list = ["JSON", "XML", "SDF", "CSV", "PNG", "TXT"]
    :return: string is used as an argument invoking request_data() function call
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    return f"{base}/{input_spec}/{operation_spec}/{output_format}"

# TODO
def build_parser():
    # parser = argparse.ArgumentParser(description="This module can help you to specify UPL parameters")
    # parser.add_argument('--async', action=argparse.BooleanOptionalAction)
    # # parser.add_argument('-d', '--domain', choices=["compound", "substance", "assay"], default="compound")
    # # parser.add_argument('-n', '--namespace', choices=["cid", "name", "smiles"], default="cid")
    # parser.add_argument('--ids', type=int, nargs='*')
    # return parser
    pass