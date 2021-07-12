import pathlib
from typing import Any, Dict, List

import click

from molecad.data.db import upload_data
from molecad.data.downloader import execute_requests
from molecad.data.utils import (
    check_dir,
    chunked,
    converter,
    file_path,
    parse_first_and_last,
    read_json,
    write_json,
)
from molecad.settings import setup


@click.group(
    help="Консольная утилита для извлечения информации о свойствах молекул с серверов Pubchem, "
    "сохранения полученных данных в файлы формата JSON, а также загрузки содержимого в базу "
    "данных MongoDB."
)
def cli():
    pass


@cli.command(
    help="Извлекает и сохраняет информацию из базы данных Pubchem – 'Compound'. "
    "Для совершения запроса к серверу необходимо уточнить диапазон идентификаторов "
    "интересуемых соединений - 'CID'. Данные извлекаются путем отправления запроса на "
    "сгенерированный URL. Из-за ограничения на количество символов в строке URL, "
    "отправка запроса, включающего все желаемые идентификаторы, может привести к ошибке, "
    "поэтому рекомендуется использовать опцию `--url-size`, которая по умолчанию равна `100`. "
    "Полученные данные сохраняются в файл, причем имеется возможность сохранять файл "
    "порциями до окончания работы программы, указав в качестве опции `--f-size` желаемое "
    "количество идентификаторов в сохраняемом файле, которое по умолчанию равна `1000`. "
)
@click.option("--start", required=True, type=int, help="Первое значение из запрашиваемых CID.")
@click.option("--stop", required=True, type=int, help="Последнее значение из запрашиваемых CID.")
@click.option(
    "--url-size",
    type=int,
    default=100,
    show_default=True,
    help="Максимальное число идентификаторов в одном запросе, по умолчанию равно 100.",
)
@click.option(
    "--f-size",
    type=int,
    default=1000,
    show_default=True,
    help="Максимальное число идентификаторов в сохраняемом файле, по умолчанию равно 1000.",
)
def fetch(start: int, stop: int, url_size: int, f_size: int) -> None:
    fetch_dir = pathlib.Path(setup.fetch_dir)
    new_dir = check_dir(fetch_dir, start, stop)
    data = execute_requests(start, stop + 1, url_size)
    chunks = chunked(data, f_size)
    for chunk in chunks:
        first, last = parse_first_and_last(chunk)
        file = file_path(new_dir, first, last)
        write_json(file, chunk)


@cli.command(
    help="Загружает содержимое файлов из указанной директории в коллекцию базы данных 'molecad' "
    "с помощью `db.collection.insert_many({..})`. MongoDB имеет ограничение на количество "
    "документов, загружаемых таким образом за один раз, которое равно 100000 документов. "
    "Если документ, который вы желаете загрузить имеет больший размер, то перед загрузкой в "
    "базу разделите его на части с помощью команды `split`"
)
@click.option(
    "--file_dir",
    required=True,
    type=pathlib.Path,
    help="Путь до директории, содержащей chunked-файлы, каждый из которых представляет собой "
    "список, длинной до 100000 элементов.",
)
@click.option(
    "--collection",
    required=True,
    type=str,
    help="Название коллекции MongoDB, в которую будут загружены файлы.",
)
def populate(f_dir: pathlib.Path, collection: str) -> None:
    n = 0
    click.echo(f"Произвожу импорт из папки {f_dir}")
    for file in f_dir.iterdir():
        click.echo(f"Импортирую файл {file}")
        data: List[Dict[str, Any]] = converter(read_json(file))
        upload_data(data, collection)
        n += len(data)
    click.echo(f"Загружено документов в коллекцию {n}")


@cli.command(
    help="При использовании команды `db.collection.insert_many({..})` имеется ограничение на "
    "максимально допустимое количество добавляемых документов за один раз равное 100000. "
    "Данная функция служит для того, чтобы разрезать JSON-файлы, превышающие указанное выше "
    "ограничение, на файлы меньшего размера. "
)
@click.option(
    "--file",
    required=True,
    type=pathlib.Path,
    help="Файл не подходящий под критерии загрузки файлов в MongoDB.",
)
@click.option(
    "--f-size",
    default=1000,
    show_default=True,
    type=int,
    help="Максимальное число идентификаторов в сохраняемом файле, по умолчанию равно 1000.",
)
def split(file: pathlib.Path, f_size: int) -> None:
    split_dir = pathlib.Path(setup.split_dir)
    data: List[Dict[str, Any]] = converter(read_json(file))
    click.echo(f"Открываю файл {file}")
    start, stop = parse_first_and_last(data)
    new_dir = check_dir(split_dir, start, stop)
    for i, chunk in enumerate(chunked(data, f_size), start=1):
        first, last = parse_first_and_last(chunk)
        ch_path = file_path(new_dir, first, last)
        write_json(ch_path, chunk)
        click.echo(f"Записываю в файл {ch_path}")


cli.add_command(fetch)
cli.add_command(populate)
cli.add_command(split)


if __name__ == "__main__":
    cli()
