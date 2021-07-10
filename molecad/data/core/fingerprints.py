from typing import Any, Dict

from mongordkit.Database import create, write, utils, registration
from mongordkit.Search import similarity, substructure, utils
from rdkit import Chem
import pymongo

from molecad.data.core.db import db, retrieve_smiles

rdkit_schema = db["test_rdkit_schema"]
mol = db["test"]


def schema_from_smiles(smiles: str) -> Dict[str, Any]:
    rdmol = Chem.MolFromSmiles(smiles)
    scheme = registration.MolDocScheme()
    return scheme.generate_mol_doc(rdmol)


def prepare_for_search():
    res = retrieve_smiles("test")
    for cid, smiles in res.items():
        schema = schema_from_smiles(smiles)
        rdkit_schema.insert_one(schema)
        mol.update_one({"CID": cid}, {"$set": {"rdkit_schema_index": schema["index"]}}, upsert=True)


def search(sub_smiles):
    q_mol = Chem.MolFromSmiles(sub_smiles)
    substructure.AddPatternFingerprints(rdkit_schema)
    index_lst = substructure.SubSearch(q_mol, rdkit_schema, chirality=False)
    # index_lst = substructure.SubSearchNaive(q_mol, rdkit_schema, chirality=False)
    cursor = mol.find({"rdkit_schema_index": {"$in": index_lst}})
    return list(cursor)


def main():
    pass


if __name__ == "__main__":
    prepare_for_search()
    res = search("CC(=O)OC")
    print([i['CanonicalSMILES'] for i in res])
