import json
import os
import uuid

DATA_DIR = os.path.expanduser("~/PycharmProjects/moleCad/data")


def save_data(obj) -> str:
    name = str(uuid.uuid4()) + ".json"
    path = os.path.join(DATA_DIR, name)
    with open(path, "wt") as f:
        json.dump(obj, f)
    return path
