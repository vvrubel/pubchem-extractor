from typing import TypeVar, Union

from molecad.types_ import (
    Domain,
    OperationAlone,
    OperationComplex,
    Out,
    PropertyTags,
)

IdT = TypeVar("IdT", int, str)
T = TypeVar("T")


def join_w_comma(*args: T) -> str:
    """
    The function formats all passed arguments from iterable type to a
    string, where ``args`` are separated by commas.
    :param args: arguments mapped as a tuple.
    :return: string of comma-separated values without whitespaces.
    """
    return ",".join(f"{i}" for i in args)


def input_specification(domain: Domain, namespace: str, identifiers: str) -> str:
    """
    This function formats specified input parameters that will be passed for
    building URL string.
    :param domain: must be a string and refers to database - compound |
    substance | assay.
    :param namespace: must be a string and can assign values depending on the
    domain. If value is complex, it passes from ``joined_namespaces()``,
    but now it is not implemented.
    :param identifiers: can be an integer or a string; the first one is
    implemented if the search is performed by single numeric id, in this
    case, the value also can be passed casted to string. In case number of
    ids is more then one (sequence), it should be previously goes to
    ``joined_identifiers()`` function.
    :return: string value formatted as "<domain>/<namespace>/<identifiers>"
    that is the first part of URL.
    """
    return f"{domain}/{namespace}/{identifiers}"


def operation_specification(
    operation: Union[OperationAlone, OperationComplex], tags: str = None
) -> str:
    """
    Dictates what to do with the input records - what information about it,
    you want to retrieve.
    The function formats the incoming parameters for further transmission to
    the URL builder.
    :param operation: describes action to execute with input identifiers.
    :param tags: string of comma-separated values without whitespaces must
    be passed, if only ``operation`` == "compound_property" or //"xrefs".
    :return: the value will be used in a ``build_url()`` function call as
    the second argument.
    """
    if tags is None:
        return f"{operation}"
    else:
        return f"{operation}/{tags}"


def build_url(input_spec: str, operation_spec: str, output_format: Out) -> str:
    """
    The command generates the URL string which is necessary to make request
    to PubChem’s PUG REST service. The URL structure consists of three
    pre-build parts that is assigned as parameters when invoking the function.
    :param input_spec: a list of identifiers for PubChem records which are
    casted to string.
    :param operation_spec: the operation to be performed with the passed
    identifiers.
    :param output_format: describes the format of the output.
    :return: string which is used as an argument of ``request_data()`` function
    call.
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    return f"{base}/{input_spec}/{operation_spec}/{output_format}"


def prepare_request(
    domain: Domain,
    namespace: str,
    identifiers: list[IdT],
    operation: Union[OperationAlone, OperationComplex],
    output: Out,
    tags: list[PropertyTags] = None,
) -> str:
    """
    Prepares all arguments to be passed to URL builder.
    Mainly, joins all list arguments with comma to string.
    :param domain: directly passes to ``build_url()`` function call.
    :param namespace: directly passes to ``build_url()`` function call.
    :param identifiers: list of string valuers or numbers which are
    associated with identifiers. It also may be passed from
    ``generate_ids()`` or  ``chunked()`` function calls.
    :param operation: if ``operation`` == "compound_property", "xrefs" or
    "target_type", the function needs additional ``tags`` argument.
    :param output: directly passes to ``build_url()`` function call.
    :param tags: must be sequence if specified.
    :return: string is used as an argument invoking ``request_data()`` function
    call.
    """
    joined_identifiers = join_w_comma(*identifiers)
    input_spec = input_specification(domain, namespace, joined_identifiers)
    if tags is None:
        operation_spec = operation_specification(operation)
    else:
        joined_tags = join_w_comma(*tags)
        operation_spec = operation_specification(operation, joined_tags)
    url = build_url(input_spec, operation_spec, output)
    return url
