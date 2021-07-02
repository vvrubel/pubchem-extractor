import json
import time
from pathlib import Path
from typing import Dict

import click
import requests
from loguru import logger

from molecad.data.db import
from molecad.data.request import (
    chunked,
    delay_iterations,
    generate_ids,
    prepare_request,
    request_data_json,
)
from molecad.types_ import (
    Domain,
    NamespCmpd,
    OperationComplex,
    Out,
    PropertyTags,
)


def execute_request(start: int, stop: int, chunk_size: int) -> Dict[int, dict]:
    """
    .. note:: В текущей версии сервиса доступен запрос свойств молекул из базы данных ``Compound``.
    Аргументы функции ``generate_ids(start, stop)`` по умолчанию равны 1 и 201 соотвественно и
    могут не указываться явно, что соответствует тестовым запросам к серверу Pubchem;
    в случае формирования базы данных эти значения должны быть явно указаны в качестве
    аргументов, как 1 и 500001 соответственно. Для последующих запросов ``start`` = ``stop`` от
    предыдущего запроса, а ``stop`` увеличивается на 500000.
    Второй аргумент, передаваемый в функцию ``chuncked`` - ``chunk_size``, в рамках для тестировых
    запросов по умолчанию равен 100 и может не указываться явно; при формировании базы данных со
    свойствами молекул должен быть равен 1000.
    """
    domain = Domain.COMPOUND
    namespace_prefix = NamespCmpd.CID
    namespace_suffix = None
    operation = OperationComplex.PROPERTY
    tags = (
        PropertyTags.MOLECULAR_FORMULA,
        PropertyTags.MOLECULAR_WEIGHT,
        PropertyTags.CANONICAL_SMILES,
        PropertyTags.INCHI,
        PropertyTags.IUPAC_NAME,
        PropertyTags.XLOGP,
        PropertyTags.H_BOND_DONOR_COUNT,
        PropertyTags.H_BOND_ACCEPTOR_COUNT,
        PropertyTags.ROTATABLE_BOND_COUNT,
        PropertyTags.ATOM_STEREO_COUNT,
        PropertyTags.BOND_STEREO_COUNT,
        PropertyTags.VOLUME_3D,
    )
    output = Out.JSON

    data = {}
    t_start = time.monotonic()
    chunks = chunked(generate_ids(start, stop), chunk_size)
    for i in delay_iterations(chunks):
        url = prepare_request(i,
                              domain,
                              namespace_prefix,
                              namespace_suffix,
                              operation,
                              tags,
                              output)
        logger.debug("Requesting URL: {}", url)
        try:
            res = request_data_json(url)
        except requests.HTTPError:
            logger.error("Error occurred: {}", exc_info=True)
            break
        else:
            logger.debug("Response content: {}", res)
            for k, v in zip(i, res['PropertyTable']['Properties']):
                data[k] = v
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    logger.info("Download took {}", t_run)
    return data

def read_data_file(inp_file: str) -> Dict[int, dict]:
    inp_file_path = Path(inp_file).resolve()
    with open(inp_file_path, 'r') as f:
        data: Dict[str, dict] = json.load(f)


@click.command()
@click.argument()
@click.option(
    '--start', default=1, prompt=True, type=int,
    help='Первое значение из запрашиваемых CID'
)
@click.option(
    '--stop', default=201, prompt=True, type=int,
    help='Последнее значение из запрашиваемых CID'
)
@click.option(
    '--chunk-size', default=100, prompt=True, type=int,
    help='Число идентификаторов в одном запросе'
)
def main():




if __name__ == '__main__':
    main()
