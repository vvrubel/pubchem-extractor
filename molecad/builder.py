from typing import Optional





def build_url(input_spec: str, operation_spec: str, output: str) -> str:
    """
    The command generates the URL to query into PubChem’s PUG REST service.
    The structure of the URL consists of three parts defined here as parameters.
    When the function invoked, pre-built argument must be specified.

    :param input_spec: string represented as "<domain>/<namespace>/<identifiers>"
    :param operation_spec: string
    :param output: string can take value from the list = ["JSON", "XML", "SDF", "CSV", "PNG", "TXT"]
    :return: string is used as an argument invoking request_data() function call
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    return f"{base}/{input_spec}/{operation_spec}/{output}"


def input_specification(domain: str, namespace: str, identifiers: str) -> str:
    """
    The function builds the first part of URL to set up one or more input identifiers representing records in one PubChem databases.
    :param domain:
    :param namespace:
    :param identifiers:
    :return: string is formatted as `<domain>/<namespace>/<identifiers>`
    """
    return f"{domain}/{namespace}/{identifiers}"


def join_w_comma(*args: object) -> str:
    """
    Format all arguments to a string.
    :param args:
    :return:
    """
    return ",".join(map(str, args))


def operation_specification(operation: str, property_tags:  Optional = None) -> str:
    if operation == "property":
        if property_tags is not None:
            return f"{operation}/{property_tags}"
        # TODO: make class MissingField
    else:
        return f"{operation}"
