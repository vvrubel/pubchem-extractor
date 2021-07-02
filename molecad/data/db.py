import json
from pathlib import Path

from loguru import logger
from pymongo import MongoClient


def data_to_db():
    client = MongoClient(URI)
    db = client[DB_NAME]
    molecad = db.test_collection
    if molecad.count_documents != 0:
        molecad.drop()
    logger.info("The collection was dropped.")

    input_dir_path = Path(INPUT_DIR).resolve()
    for i, file_path in enumerate(input_dir_path.iterdir()):
        logger.info(f"Loading data from {file_path}. {i * 1000} are done.")
        with open(file_path, "rt") as f:
            data = json.load(f)
        molecad.insert_many(data)


if __name__ == '__main__':
    main()