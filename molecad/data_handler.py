import click
import json
import uuid

from loguru import logger
from pathlib import Path
from typing import Any, Generator


def check_dir(out_dir: Path) -> None:
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except (FileNotFoundError, FileExistsError) as e:
        logger.error(f"Error occurred: {e}", exc_info=True)
        raise


def save_data_json(obj: Any, out_dir: Path) -> Path:
    """
    Saves json file outputted from ``request_data_json`` function call.
    :param obj: JSON object.
    :param out_dir: output directory.
    :return: path to json file.
    """
    name = str(uuid.uuid4()) + '.json'
    f_path = Path(out_dir) / name
    with open(f_path, 'wt') as f:
        json.dump(obj, f, indent=2)
    return f_path


def cutter(inp_file: Path, out_dir: Path, maxsize: int = 1000) -> Generator[Path, None, None]:
    from pubchem_requests import chunked
    with open(inp_file, 'r') as f:
        data: dict[str, dict] = json.load(f)
    check_dir(out_dir)
    for i, chunk in enumerate(chunked(data.values(), maxsize=maxsize)):
        chunk_path = save_data_json(chunk, out_dir)
        logger.info(f"Data for {i + 1} is saved in {chunk_path}")
        yield chunk_path


@click.command()
@click.argument('inp_file')
@click.argument('out_dir')
@click.option(
    '--maxsize', '-m',
    help="how many documents should be in a chunked json"
)
def main(inp_file: Path, out_dir: Path, maxsize: int = 1000):
    inp_file = Path(inp_file).resolve()
    out_dir = Path(out_dir).resolve()
    file_list = []
    for i in cutter(inp_file, out_dir):
        file_list.append(i)
    logger.info(f"Data saved in {out_dir}.")


if __name__ == '__main__':
    main()
