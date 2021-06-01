import requests

from .types_ import *

BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
OUTPUT_FORMAT = "JSON"
# OPTIONAL = "{[?<operation_options>]}"


def build_url(input_spec: str, operation_spec: str) -> str:
    return f"{BASE_URL}/{input_spec}/{operation_spec}/{OUTPUT_FORMAT}"


def input_specification(domain, namespace, identifiers):
    return f"{domain}/{namespace}/{identifiers}"


def build_identifiers(*identifiers):
    return ",".join(map(str, identifiers))


def operation_specification():
    '''
    If domain == compaund
    :return:
    '''


def request_data(url, **params):
    response = requests.get(url, params=params)
    data = response.json()
    return data
