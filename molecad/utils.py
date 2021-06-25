import click
import json
import uuid

from loguru import logger
from pathlib import Path
from typing import Any, Generator

from molecad.downloader import chunked


def save_data_json(obj: Any) -> Path:
    """
    Saves outputted json file.
    :param obj: JSON object.
    :param num: prefix for naming.
    :return: path to json file.
    """
    try:
        d_path = Path().resolve() / 'data'
        d_path.mkdir(parents=True, exist_ok=True)
    except (FileNotFoundError, FileExistsError) as e:
        logger.error(f"Error occurred: {e}", exc_info=True)
    else:
        name = str(uuid.uuid4()) + '.json'
        f_path = Path(d_path) / name
        with open(f_path, 'wt') as f:
            json.dump(obj, f, indent=2)
    return f_path


def cutter(big_path: Path, maxsize: int = 1000) -> Generator[Path, None, None]:
    with open(big_path, 'r') as f:
        big_json = json.load(f)
    records = list(big_json.values())
    for chunk in chunked(records, maxsize=maxsize):
        chunk_path = save_data_json(chunk)
        logger.info(f"Data for {chunk} is saved in {chunk_path}")
        yield chunk_path


@click.command()
@click.argument('big')
@click.option(
    '--maxsize', '-m',
    help="how many documents should be in a chunked json"
)
def main(big: Path, maxsize: int = 1000):
    big_path = Path(big).resolve()
    file_list = []
    for i in cutter(big_path):
        file_list.append(i)
        logger.info(f"Data saved in {i}.")


if __name__ == '__main__':
    main()
