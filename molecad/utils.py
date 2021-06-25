import json
from loguru import logger
from pathlib import Path
from typing import Generator

from molecad.downloader import save_data_json


def cutter(f_path: Path, maxsize: int = 1000) -> Generator[Path, None, None]:
    with open(f_path, "r") as f:
        big_json = json.load(f)
    records = list(big_json.values())
    for chunk in chunked(records, maxsize=maxsize):
        chunk
        yield chunk_path



json_file = Path('/data/c47b4dbd-68a2-4670-bdac-9f9d0f983b1d.json').resolve()
for file_name in cutter(json_file):
    logger.info("Data saved in {}", file_name)


if __name__ == '__main__':
    main()