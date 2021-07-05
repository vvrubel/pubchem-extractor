import pathlib
from typing import Any, Dict, List

import click

from molecad.data.core.db import drop_collection, upload_data
from molecad.data.core.downloader import execute_request
from molecad.data.core.utils import check_dir, chunked, file_name, read, write


@click.group()
def cli():
    pass


@click.command()
@click.option(
    "--out-dir",
    required=True,
    type=pathlib.Path,
    help="Путь до директории, в которую будет записан файл",
)
@click.option(
    "--start", default=1, required=True, type=int, help="Первое значение из запрашиваемых CID"
)
@click.option(
    "--stop", default=201, required=True, type=int, help="Последнее значение из запрашиваемых CID"
)
@click.option(
    "--size", default=100, required=True, type=int, help="Число идентификаторов в одном запросе"
)
def fetch(out_dir: pathlib.Path, start: int, stop: int, size: int) -> None:
    """
    Выполняет запрос и пишет данные в файл.
    :param out_dir: Папка, в которую будет записан файл.
    :param start: Первое значение из запрашиваемых CID.
    :param stop: Последнее значение из запрашиваемых CID.
    :param size: Число идентификаторов в одном запросе.
    """
    data = execute_request(start, stop, size)
    check_dir(out_dir)
    file = file_name(out_dir)
    write(file, data)


@click.command()
@click.option("--file", required=True, type=pathlib.Path)
@click.option("--f-dir", required=True, type=pathlib.Path)
@click.option("--size", default=1, type=int)
def split(file: pathlib.Path, f_dir: pathlib.Path, size: int) -> None:
    """
    Разрезает большой JSON на чанки.
    :param file: путь до большого JSON.
    :param f_dir: путь до директории, в которую будут записаны созданные файлы.
    :param size: размер чанка.
    """
    check_dir(f_dir)
    data: List[Dict[str, Any]] = list(read(file).values())
    click.echo(f"Открываю файл {file}")
    for i, chunk in enumerate(chunked(data, size), start=1):
        ch_path = file_name(f_dir)
        write(ch_path, chunk)
        click.echo(f"Записываю в файл {ch_path}")


@click.command()
@click.option("--f-dir", required=True, type=pathlib.Path)
def populate(f_dir: pathlib.Path):
    n = 0
    drop_collection()
    click.echo(f"Произвожу импорт из папки {f_dir}")
    for file in f_dir.iterdir():
        click.echo(f"Импортирую файл {file}")
        data = read(file)
        upload_data(data)
        n += len(data)
    click.echo(f"Загружено документов в коллекцию {n}")


cli.add_command(fetch)
cli.add_command(split)
cli.add_command(populate)


if __name__ == "__main__":
    cli()
