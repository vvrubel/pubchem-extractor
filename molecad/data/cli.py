import json
import sys
import time
from pathlib import Path
from typing import Dict, Iterator

import click
import click_completion
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
from molecad.errors import (
    DirNotFoundError
)
from molecad.types_ import (
    DataHandler,
    Domain,
    NamespCmpd,
    OperationComplex,
    Out,
    PropertyTags,
)

click_completion.init()


def check_dir(out_dir: str) -> None:
    """
    Проверяет существует ли директория, и создает ее в случае отсутствия.
    :param out_dir: строкоевое представления пути до директории.
    :return: None
    """
    out_dir_path = Path(out_dir).resolve()
    try:
        out_dir_path.mkdir(parents=True, exist_ok=False)
    except DirNotFoundError as e:
        logger.exception("Директория уже существует")


def cutter(obj, maxsize: int = 100) -> Iterator[Path]:
    from molecad.data.request import chunked
    for i, chunk in enumerate(chunked(obj.values(), maxsize=maxsize)):
        chunk_path = save_data_json(chunk)
        logger.info(f"Data for {i + 1} is saved in {chunk_path}")
        yield chunk


@click.group(chain=True, invoke_without_command=True)
@click.argument('action', help='Способ которым будут обработаны данные')
@click.pass_context
def cli(ctx: click.Context, action: DataHa):
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        sys.exit(0)
    else:
        if ctx.invoked_subcommand == "execute_request":
            click.echo("Начинаю выполнять запросы к базам данных Pubchem")
            ctx.obj = execute_request()
            if action == DataHandler.SAVE:



@click.command()
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

    results = {}
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
                results[k] = v
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    logger.info("Download took {}", t_run)
    return results


@click.command()
@click.argument('inp_file', prompt=True, type=click.Path(exists=True))
def read_data_file(inp_file: str) -> Dict[int, dict]:
    inp_file_path = Path(inp_file).resolve()
    with open(inp_file_path, 'r') as f:
        data: Dict[str, dict] = json.load(f)



@cli.command()
@click.argument(
    'out_dir', default='./molecad/data/response_files', prompt=True,
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    help='Путь до директории, в которуюбудет сохранен файл'
)
@click.option(

)
 @click.pass_obj
def save_data_json(obj: Dict[int, dict], out_dir: str) -> Path:
    """
    Сохраняет данные, полученные от сервера путем выолнения функции ``request_data_json``, в файл.
    :param obj: JSON объект.
    :param out_dir: путь до директории, в которую будут сохранены файлы, указанный относительно
    папки проекта, например: ``./molecad/data/response_files``
    :return: путь до сохраненного файла.
    """
    check_dir(out_dir)
    file_list = []
    for i in cutter(obj):
        file_list.append(i)
    logger.info(f"Data saved in {out_dir}."
    name = str(obj.keys()) + '.json'
    f_path = Path(out_dir) / name
    with open(f_path, 'wt') as f:
        json.dump(obj.values(), f, indent=2)
    return f_path


# @click.command()
# @click.argument('inp_file')
# @click.option(
#     '--out-dir', show_envvar=True, default='./molecad/data/response_files', prompt=True,
#     type=click.Path(exists=True, dir_okay=True, resolve_path=True),
#     help='Путь до директории, в которуюбудет сохранен файл'
# )
# @click.option(
#     '--chunk-size', default=100, prompt=True, type=int,
#     help='Число идентификаторов в одном запросе'
# )
# def main(inp_file: str, out_dir: str, maxsize: int = 1000):

cli.add_command(execute_request)
cli.add_command(read_data_file)
cli.add_command(save_data_json)


if __name__ == '__main__':
    cli()