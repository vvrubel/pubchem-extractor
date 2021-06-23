from typing import Iterable, Optional, Sequence, Union

from molecad.types_ import (
    IdT,
    InputDomains,
    InputNamespaces,
    MissingValue,
    Operations,
    OutputFormats,
    PropertyTags,
    SearchPrefix,
    SearchSuffix,
    T
)


def join_w_comma(*args: IdT) -> str:
    """
    Formats all passed arguments from iterable to a string separated by commas
    :param args: arguments mapped as a tuple
    :return: string of comma-separated values without whitespaces.
    """
    return ",".join(f'{i}' for i in args)


def joined_namespaces(
        prefix: Union[str, SearchPrefix],
        suffix: Optional[SearchSuffix]
) -> str:
    """
    Prepares namespace parts which are needed to be concatenated with "/".
    :param prefix: value from ``structure search`` or  ``fast search``
    namespace subclass.
    :param suffix: if value is None, error must be raised.
    :return: formatted a namespace is passed to ``input_specification``
    function call.
    """
    if suffix is not None:
        return f"{prefix}/{suffix}"
    else:
        raise MissingValue


def input_specification(
        domain: Union[str, InputDomains],
        namespace: Union[str, InputNamespaces],
        identifiers: str,
) -> str:
    """
    This function will be used to format specified input parameters that
    will be passed for building URL string.
    :param domain: string value refers to database - compound | substance |
    assay.
    :param namespace: must be string and can assign values depending on the
    domain. If value is complex, it passed from ``joined_namespaces``.
    :param identifiers: passes from ``joined_identifiers``.
    :return: string value that  is the first part of URL and formatted
    as "<domain>/<namespace>/<identifiers>"
    """
    return f"{domain}/{namespace}/{identifiers}"


def operation_specification(
        operation: str,
        tags: Optional[str] = None
) -> str:
    """
    The function dictates what to do with the input records - what
    information about them you want to retrieve. The function will be used
    to format specified operation parameters that will be passed for
    building URL string.
    :param operation: describes action to execute with input identifiers.
    :param tags: string of comma-separated values without whitespaces must
    be passed, if only ``operation`` == "compound_property", "xrefs" or
    "target_type".
    :return: the value will be used in a ``build_url()`` function call as 
    the second argument.
    """
    if tags is None:
        return f"{operation}"
    else:
        return f"{operation}/{tags}"


def build_url(
        input_spec: str,
        operation_spec: str,
        output_format: str
) -> str:
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
        domain: InputDomains,
        namespace: Union[str, InputNamespaces],
        identifiers: list[IdT],
        operation: Operations,
        output: OutputFormats,
        tags: Optional[Sequence[PropertyTags]] = None
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
