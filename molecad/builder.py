from typing import Optional, Sequence

from molecad.types_ import (
    # AssayInputNamespaces,
    # CompoundInputNamespaces,
    # FastSearchPrefix,
    # GeneralInputNamespaces,
    IdT,
    InputDomains,
    # NamespaceSearchSuffix,
    Operations,
    OutputFormats,
    PropertyTags,
    # StructureSearchPrefix,
    # SubstanceInputNamespaces,
    # TargetTypes,
    # Xrefs
)


def join_w_comma(*args: object) -> str:
    return ",".join(f'{i}' for i in args)


def namespace_concatinator(
        prefix: str,
        suffix: Optional[str] = None
) -> str:
    if suffix is None:
        return f"{prefix}"
    else:
        return f"{prefix}/{suffix}"


def input_specification(
        domain: str,
        namespace: str,
        identifiers: str
) -> str:
    """
    This function will be used to format specified input parameters that
    will be passed for building URL string.
    :param domain:
    :param namespace:
    :param identifiers:
    :return:
    """
    return f"{domain}/{namespace}/{identifiers}"


def operation_specification(
        operation: str,
        tags: Optional[str] = None
) -> str:
    """
    This function will be used to format specified operation parameters that
    will be passed for building URL string.
    :param operation:
    :param tags:
    :return:
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
    :param input_spec:
    :param operation_spec:
    :param output_format:
    :return:
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    return f"{base}/{input_spec}/{operation_spec}/{output_format}"


def prepare_request(
        domain: InputDomains,
        namespace: str,
        identifiers: IdT,
        operation: Operations,
        output: OutputFormats,
        tags: Optional[Sequence[PropertyTags]] = None
) -> str:
    joined_identifiers = join_w_comma(*identifiers)
    input_spec = input_specification(domain, namespace, joined_identifiers)
    if tags is None:
        operation_spec = operation_specification(operation)
    else:
        joined_tags = join_w_comma(*tags)
        operation_spec = operation_specification(operation, joined_tags)
    url = build_url(input_spec, operation_spec, output)
    return url
