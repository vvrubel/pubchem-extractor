import pathlib
from typing import Any, Dict, List

import click

from molecad.data.db import drop_collection, upload_data
from molecad.data.downloader import execute_requests
from molecad.data.utils import check_dir, chunked, converter, file_name, read, write
from molecad.settings import setup


@click.group(
    help="Утилита для извлечения информации с серверов Pubchem и ее записи в файлы или в "
    "локальную базу данных с помощью интерфейса командной строки."
)
def cli():
    pass


@click.command(
    help="Выполняет запрос к серверу Pubchem, извлекает данные из ответа сервера и пишет их в файл."
)
@click.option(
    "--start", required=True, type=int, help="Первое значение из запрашиваемых CID."
)
@click.option(
    "--stop", required=True, type=int, help="Последнее значение из запрашиваемых CID."
)
@click.option(
    "--req-size",
    default=100,
    type=int,
    help="Максимальное число идентификаторов в одном запросе.",
)
@click.option(
    "--f-size",
    default=1000,
    type=int,
    help="Максимальное число идентификаторов в сохраняемом файле.",
)
def fetch(start: int, stop: int, req_size: int, f_size: int) -> None:
    fetch_dir = pathlib.Path(setup.fetch_dir)
    new_dir = check_dir(fetch_dir, start, stop)
    data = execute_requests(start, stop, req_size)
    chunks = chunked(data, f_size)
    for chunk in chunks:
        file = file_name(new_dir)
        write(file, chunk)


@click.command(
    help="Разрезает большой JSON на чанки меньшего размера для последующей загрузки в "
    "MongoDB, что необходимо из-за внутренних ограничений MongoDB на количество "
    "документов, загружаемых за один раз одним файлом."
)
@click.argument("file", type=pathlib.Path)
@click.option(
    "--f-dir",
    required=True,
    type=pathlib.Path,
    help="Путь до директории, в которую будут записаны созданные chunked-файлы – не должна "
    "существовать до начала выполнения записи файлов.",
)
@click.option(
    "--size", default=1000, type=int, help="Максимальное число элементов в одном chunked-файле"
)
def split(file: pathlib.Path, f_dir: pathlib.Path, size: int) -> None:
    check_dir(f_dir)
    data: List[Dict[str, Any]] = converter(read(file))
    click.echo(f"Открываю файл {file}")
    for i, chunk in enumerate(chunked(data, size), start=1):
        ch_path = file_name(f_dir)
        write(ch_path, chunk)
        click.echo(f"Записываю в файл {ch_path}")


@click.command(help="Загружает chunked-файлы из указанной директории в локальную базу MongoDB.")
@click.option(
    "--f-dir",
    required=True,
    type=pathlib.Path,
    help="Путь до директории, содержащей chunked-файлы, содержимое каждого из которых "
    "представляет собой список, длинной до 100000 элементов - ограничение MongoDB",
)
@click.option(
    "--collection",
    required=True,
    type=str,
    help="Название коллекции MongoDB, в которую будут загружены файлы.",
)
def populate(f_dir: pathlib.Path, collection: str) -> None:
    n = 0
    drop_collection(collection)
    click.echo(f"Произвожу импорт из папки {f_dir}")
    for file in f_dir.iterdir():
        click.echo(f"Импортирую файл {file}")
        data: List[Dict[str, Any]] = converter(read(file))
        upload_data(data, collection)
        n += len(data)
    click.echo(f"Загружено документов в коллекцию {n}")


cli.add_command(fetch)
cli.add_command(split)
cli.add_command(populate)


if __name__ == "__main__":
    cli()
