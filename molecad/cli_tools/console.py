import pathlib
from typing import Any, Dict, List

import click
import pymongo
import pymongo.errors
import rdkit.RDLogger
from loguru import logger
from mongordkit import Search

from molecad.cli_tools.db import (
    create_indexes,
    create_molecule,
    delete_broken,
    drop_db,
    upload_data,
)
from molecad.cli_tools.downloader import execute_requests
from molecad.settings import Settings, settings
from molecad.utils import (
    chunked,
    converter,
    create_dir,
    file_path,
    parse_first_and_last,
    read_json,
    timer,
    write_json,
)

rdkit.RDLogger.DisableLog("rdApp.*")


@click.group(
    help="Консольная утилита для извлечения информации о свойствах молекул с серверов Pubchem, "
    "сохранения полученной информации в файлы формата JSON. Также с помощью утилиты можно "
    "загружать полученные данные в базу MongoDB и подготавливать их для поиска молекул. ",
    invoke_without_command=True,
)
@click.option(
    "--env",
    type=click.Choice(["prod", "dev"], case_sensitive=False),
    help="Опция позволяет выбирать конфигурационный файл, содержащий переменные окружения",
)
@click.pass_context
def molecad(ctx: click.Context, env: str):
    if env == "prod":
        ctx.obj = Settings(_env_file=".env.prod", _env_file_encoding="utf-8")
    else:
        ctx.obj = settings

    settings.setup_logging()
    logger.trace(f"Команда запущена с контекстом {ctx.obj}.")


@molecad.command(
    help="Извлекает и сохраняет информацию из базы данных Pubchem – 'Compound'. "
    "Для совершения запроса к серверу необходимо уточнить диапазон идентификаторов "
    "интересуемых соединений - 'CID'. Данные извлекаются путем отправления запроса на "
    "сгенерированный URL. Из-за ограничения на количество символов в строке URL, "
    "отправка запроса, включающего все желаемые идентификаторы, может привести к ошибке, "
    "поэтому рекомендуется использовать опцию `--size`, которая по умолчанию равна `100`. "
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
    type=int,
    help="Максимальное число идентификаторов в одном запросе, по умолчанию равно 100.",
)
@click.option(
    "--f-size",
    default=1000,
    show_default=True,
    type=int,
    help="Максимальное число идентификаторов в сохраняемом файле, по умолчанию равно 1000.",
)
@click.pass_obj
@timer
def fetch(obj: Settings, start: int, stop: int, size: int, f_size: int) -> None:
    fetch_dir = pathlib.Path(obj.fetch_dir)
    new_dir: pathlib.Path = create_dir(fetch_dir, start, stop)
    logger.info(f"Файлы будут сохранены в директорию {new_dir}")

    data = execute_requests(start, stop + 1, size)
    chunks = chunked(data, f_size)
    for chunk in chunks:
        first, last = parse_first_and_last(chunk)
        file: pathlib.Path = file_path(new_dir, first, last)
        write_json(file, chunk)
        logger.success(f"Данные для CID: {first} – {last} сохранены")


@molecad.command(
    help='При использовании команды MongoDB "db.collection.insert_many({...})" имеется '
    "ограничение на максимально допустимое количество добавляемых документов за раз – 100000. "
    "Данная функция служит для того, чтобы разрезать JSON-файлы, превышающие указанное "
    "ограничение, на файлы меньшего размера с целью сделать их пригодными для загрузки в базу."
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
    help="Максимальное число идентификаторов в каждом из сохраняемых файлов, по умолчанию – 1000.",
)
@click.pass_obj
@timer
def split(obj: Settings, file: pathlib.Path, f_size: int) -> None:
    split_dir = pathlib.Path(obj.split_dir)
    data: List[Dict[str, Any]] = converter(read_json(file))
    logger.info(f"Для разделения получен файл {file.name}")

    start, stop = parse_first_and_last(data)
    new_dir: pathlib.Path = create_dir(split_dir, start, stop)
    logger.info(f"Файлы будут сохранены в директорию {new_dir}")

    for i, chunk in enumerate(chunked(data, f_size), start=1):
        first, last = parse_first_and_last(chunk)
        ch_path: pathlib.Path = file_path(new_dir, first, last)
        write_json(ch_path, chunk)
        logger.success(f"Сохранен файл {ch_path.name}")


@molecad.command(
    help="Файлы из указанной директории содержат данные скачанные из Pubchem. Количество "
    "документов в файле не должно превышать 100000, иначе данный файл будет пропущен при загрузке. "
    "В ходе работы программы для каждого загружаемого в коллекцию `properties` документа "
    "генерируется схема молекулы в коллекции `molecules`, после чего обе коллекции "
    "подготавливаются для подструкрного поиска. При указании опции `--drop` все коллекции будут "
    "удалены и созданы заново, после чего на них будет установлен уникальный индекс `CID` и "
    "индекс `index`."
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
    help='Если опция "--drop" указана, то база будет очищена.',
)
@click.pass_obj
@timer
def populate(obj: Settings, f_dir: pathlib.Path, drop: bool) -> None:
    db = pymongo.MongoClient(obj.mongo_uri)[obj.db_name]
    if drop:
        drop_db(db)

    create_indexes(db[obj.properties], db[obj.molecules])

    total, succeed, failed, created = (0 for _ in range(4))
    logger.info(f"Импорт будет производится директории {f_dir.name} в базу {db.name}")
    for file in f_dir.glob("./**/*.json"):
        logger.debug(f"Импортирую файл {file.name}")
        try:
            data: List[Dict[str, Any]] = converter(read_json(file))
            total += len(data)

            processed, inc_created = create_molecule(data, db[obj.molecules])
            inc_succeed, inc_failed = upload_data(processed, db[obj.properties])

            created += inc_created
            succeed += inc_succeed
            failed += inc_failed
        except pymongo.errors.BulkWriteError as e:
            click.secho(e, fg="red")
            continue
        except pymongo.errors.DuplicateKeyError as e:
            click.secho(e, fg="red")
            continue
    without_smiles, without_scheme = delete_broken(db[obj.properties])

    logger.success(
        f"Выбрана база данных: {db.name}\n"
        f"Коллекция представлений молекул: {db[obj.molecules].name}\n"
        f"Создано схем: {created}\n"
        f"Коллекция свойств молекул: {db[obj.properties].name}\n"
        f"Общее количество обработанных документов: {total}\n"
        f"Новых документов: {succeed},\n"
        f"Пропущено дубликатов: {failed},\n"
        f"Удалено документов без SMILES: {without_smiles},\n"
        f"Удалено документов без схемы: {without_scheme}",
    )

    logger.info("Начинаю генерировать Fingerprints")
    Search.AddPatternFingerprints(db[obj.molecules])
    logger.success("Команда AddPatternFingerprints выполнена.")
    Search.AddMorganFingerprints(db[obj.molecules], db[obj.mfp_counts])
    logger.success("Команда AddMorganFingerprints выполнена.")


molecad.add_command(fetch)
molecad.add_command(split)
molecad.add_command(populate)


if __name__ == "__main__":
    molecad()
