import json
import requests
import time
import uuid
from loguru import logger
from pathlib import Path
# from typing import Sequence

from molecad.builder import (
    prepare_request
)
from molecad.sample_sync_requests import (
    chunked,
    delay_iterations,
    generate_ids,
    request_data_json,
)
from molecad.types_ import (
    InputDomains,
    CompoundInputNamespaces,
    Operations,
    PropertyTags,
    OutputFormats,
)


def execute_request(url: str, params: dict[str, str]) -> dict[str, object]:
    logger.info("executing request {}", url)
    data = request_data_json(url, **params)
    return data


def save_data_json(obj: dict[str, object]) -> Path:
    d_path = Path().resolve() / "downloaded"
    if not Path(d_path).exists():
        d_path.mkdir()
    name = str(uuid.uuid4()) + ".json"
    f_path = Path(d_path) / name
    with open(f_path, "wt") as f:
        json.dump(obj, f, indent=2)
    return f_path


def main():
    domain = InputDomains.COMPOUND
    namespace = CompoundInputNamespaces.CID
    operation = Operations.PROPERTY
    tags = (
        PropertyTags.MOLECULAR_FORMULA,
        PropertyTags.MOLECULAR_WEIGHT,
        PropertyTags.CANONICAL_SMILES,
        PropertyTags.INCHI,
        PropertyTags.IUPAC_NAME,
        PropertyTags.XLOGP,
        PropertyTags.H_BOND_DONOR_COUNT,
        PropertyTags.H_BOND_ACCEPTOR_COUNT,
        PropertyTags.ROTATABLE_BOND_COUNT,
        PropertyTags.ATOM_STEREO_COUNT,
        PropertyTags.BOND_STEREO_COUNT,
        PropertyTags.VOLUME_3D,
    )
    output = OutputFormats.JSON

    results = {}
    t_start = time.monotonic()
    for i in delay_iterations(chunked(generate_ids(), 100), 60.0, 400):
        url = prepare_request(domain, namespace, i, operation, output, tags)
        logger.debug("Requesting URL: {}", url)
        try:
            res = execute_request(url, {})
        except requests.HTTPError:
            logger.error("Error occurred: {}", exc_info=True)
            break
        else:
            logger.debug("Response content: {}", res)
            for k, v in zip(i, res['PropertyTable']['Properties']):
                results[k] = v
    t_stop = time.monotonic()
    t_run = t_stop - t_start
    logger.info("Download took {}", t_run)

    res_file = save_data_json(results)
    logger.info("Data saved in {}", res_file)


if __name__ == '__main__':
    main()
