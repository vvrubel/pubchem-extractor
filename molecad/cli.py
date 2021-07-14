import pathlib
from typing import Any, Dict, List

import click
import pymongo.errors

from .db import create_index, delete_without_smiles, upload_data
from .downloader import execute_requests
from .errors import CreateIndexError
from .settings import Settings
from .utils import (
    check_dir,
    chunked,
    converter,
    file_path,
    parse_first_and_last,
    read_json,
    write_json,
)


@click.group(
    help="Консольная утилита для извлечения информации о свойствах молекул с серверов Pubchem, "
         "сохранения полученной информации в файлы формата JSON. Также с помощью утилиты можно "
         "загружать полученные данные в базу MongoDB и подготавливать их для поиска молекул."
         "Для выполнения команд, которые устанавливают соединение с базой данных, необходимо "
         "указать опцию --prod/--dev для выбора рабочей базы.",
    invoke_without_command=True
)
@click.option(
    '--setup',
    type=click.Choice(["PROD", "DEV"], case_sensitive=False),
    help="Опция, позволяющая выбирать конфигурационный файл, содержащий переменные окружения, "
         "который также определяет настройки базы данных."
)
@click.pass_context
def molecad(ctx: click.Context, setup: str):
    if setup == "DEV":
        ctx.obj = Settings()
    else:
        ctx.obj = Settings(_env_file='prod.env', _env_file_encoding='utf-8')
    click.secho(f"Выбрана база {ctx.obj.db_name}", fg="blue")


@molecad.command(
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
    "--size",
    default=100,
    show_default=True,
    help="Максимальное число идентификаторов в одном запросе, по умолчанию равно 100.",
)
@click.option(
    "--f-size",
    default=1000,
    show_default=True,
    help="Максимальное число идентификаторов в сохраняемом файле, по умолчанию равно 1000.",
)
@click.pass_obj
def fetch(obj: Any, start: int, stop: int, size: int, f_size: int) -> None:
    fetch_dir = pathlib.Path(obj.fetch_dir)
    new_dir = check_dir(fetch_dir, start, stop)
    click.echo(f"Файлы будут сохранены в директорию {new_dir}")
    data = execute_requests(start, stop + 1, size)
    chunks = chunked(data, f_size)
    for chunk in chunks:
        first, last = parse_first_and_last(chunk)
        click.echo(f"Сохраняю данные для CID: {first} – {last}")
        file = file_path(new_dir, first, last)
        write_json(file, chunk)


@molecad.command(
    help="При использовании команды `db.collection.insert_many({..})` имеется ограничение на "
    "максимально допустимое количество добавляемых документов за один раз равное 100000. "
    "Данная функция служит для того, чтобы разрезать JSON-файлы, превышающие указанное выше "
    "ограничение, на файлы меньшего размера."
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
    help="Максимальное число идентификаторов в сохраняемом файле, по умолчанию равно 1000.",
)
@click.pass_obj
def split(obj: Any, file: pathlib.Path) -> None:
    split_dir = pathlib.Path(obj.split_dir)
    data: List[Dict[str, Any]] = converter(read_json(file))
    click.echo(f"Открываю файл {file}")
    start, stop = parse_first_and_last(data)
    new_dir = check_dir(split_dir, start, stop)
    for i, chunk in enumerate(chunked(data, obj.f_size), start=1):
        first, last = parse_first_and_last(chunk)
        ch_path = file_path(new_dir, first, last)
        write_json(ch_path, chunk)
        click.echo(f"Записываю в файл {ch_path}")


@molecad.command(
    help="Из указанной директории загружает файлы в коллекцию MongoDB, при условии что файл "
         "содержит не более 100000 документов. Перед загрузкой документов на коллекции будут "
         "созданы уникальные индекс 'CID'."
)
@click.option(
    "--f-dir",
    required=True,
    type=pathlib.Path,
    help="Путь до директории, содержащей chunked-файлы, каждый из которых представляет собой "
         "список длиной до 100000 элементов.",
)
@click.option(
    "--drop",
    is_flag=True,
    help="Если опция '--drop' указана, то коллекция будет очищена.",
)
@click.pass_obj
def populate(obj: Any, f_dir: pathlib.Path, drop: bool) -> None:
    db = obj.get_db
    collection = db[obj.pubchem]
    if drop:
        collection.drop()
        click.secho(f"Коллекция {collection.name} была очищена.", fg='yellow')

    try:
        create_index(collection, "CID")
        click.secho(f"На коллекции {collection.name} создан индекс – \"CID\".", fg='green')
    except CreateIndexError:
        click.secho(f"На коллекции {collection.name} индекс \"CID\" уже создан.", fg='red')

    succeed = 0
    failed = 0
    click.secho(f"Произвожу импорт из папки {f_dir}", fg='blue')
    for file in f_dir.iterdir():
        try:
            click.echo(f"Импортирую файл {file}")
            data: List[Dict[str, Any]] = converter(read_json(file))
            s, f = upload_data(data, collection)
            succeed += s
            failed += f
        except pymongo.errors.BulkWriteError as e:
            click.echo(e)

    deleted = delete_without_smiles(collection)
    click.secho(f"Загружено: {succeed},\nПропущено: {failed},\nУдалено: {deleted}", fg='blue')


molecad.add_command(fetch)
molecad.add_command(populate)
molecad.add_command(split)


if __name__ == "__main__":
    molecad()
