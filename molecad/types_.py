from enum import Enum
from loguru import logger
from typing import TypeVar


IdT = TypeVar("IdT", int, str)
T = TypeVar("T")


class WrongValue(ValueError, TypeError):
    def __init__(self, message: str = "") -> None:
        super().__init__(message)
        self.message = message
        logger.error("Error occurred: {}. Value is wrong.", exc_info=True)


class Domain(Enum):
    COMPOUND = "compound"
    SUBSTANCE = "substance"
    ASSAY = "assay"


class NamespaceComp(Enum):
    """
    The class describes the namespaces that must be specified in input
    specification. The value of class elements is a tuple containing the
    type of search parameter (the first) and its group of the complexity (the
    second). The second dictates if namespace has a suffix or not.
    If ``group`` == 0, the result namespace value is a single word, else it
    must has a suffix chosen from ``Xref`` or ``SearchSuffix`` classes.
    """
    CID = "cid", 0
    NAME = "name", 0
    SMILES = "smiles", 0
    SDF = "sdf", 0
    INCHI = "inchi", 0
    INCHIKEY = "inchikey", 0
    LISTKEY = "listkey", 0
    FORMULA = "formula", 0
    FAST_FORMULA = "fastformula", 0
    SUBSTRUCTURE = "substructure", 1
    SUPERSTRUCTURE = "superstructure", 1
    SIMILARITY = "similarity", 1
    IDENTITY = "identity", 1
    FAST_IDENTITY = "fastidentity", 1
    FAST_SIMILARITY_2D = "fastsimilarity_2d", 1
    FAST_SIMILARITY_3D = "fastsimilarity_3d", 1
    FAST_SUBSTRUCTURE = "fastsubstructure", 1
    FAST_SUPERSTRUCTURE = "fastsuperstructure", 1
    XREF = "xref", 2

    def __init__(self, search: str, group: int) -> None:
        self._search_ = search
        self._group_ = group

    @property
    def search(self):
        return self._search_

    @property
    def group(self):
        return self._group_


class SearchSuffix(Enum):
    SMILES = "smiles"
    INCHI = "inchi"
    SDF = "sdf"
    CID = "cid"


class Xref(Enum):
    REGISTRY_ID = "RegistryID"
    RN = "RN"
    PUBMED_ID = "PubMedID"
    MMDB_ID = "MMDBID"
    PROTEIN_GI = "ProteinGI"
    NUCLEOTIDE_GI = "NucleotideGI"
    TAXONOMY_ID = "TaxonomyID"
    MIM_ID = "MIMID"
    GENE_ID = "GeneID"
    PROBE_ID = "ProbeID"
    PATENT_ID = "PatentID"


class Operations(Enum):
    RECORD = "record", 0
    SYNONYMS = "synonyms", 0
    SIDS = "sids", 0
    CIDS = "cids", 0
    AIDS = "aids", 0
    ASSAY_SUMMARY = "assaysummary", 0
    CLASSIFICATION = "classification", 0
    DESCRIPTION = "description", 0
    CONFORMERS = "conformers", 0
    PROPERTY = "property", 1
    XREFS = "xrefs", 2

    def __init__(self, action: str, group: int) -> None:
        self._action_ = action
        self._group_ = group

    @property
    def action(self):
        return self._action_

    @property
    def group(self):
        return self._group_


class PropertyTags(Enum):
    MOLECULAR_FORMULA = "MolecularFormula"
    MOLECULAR_WEIGHT = "MolecularWeight"
    CANONICAL_SMILES = "CanonicalSMILES"
    ISOMERIC_SMILES = "IsomericSMILES"
    INCHI = "InChI"
    INCHI_KEY = "InChIKey"
    IUPAC_NAME = "IUPACName"
    TITLE = "Title"
    XLOGP = "XLogP"
    EXACT_MASS = "ExactMass"
    MONOISOTOPIC_MASS = "MonoisotopicMass"
    TPSA = "TPSA"
    COMPLEXITY = "Complexity"
    CHARGE = "Charge"
    H_BOND_DONOR_COUNT = "HBondDonorCount"
    H_BOND_ACCEPTOR_COUNT = "HBondAcceptorCount"
    ROTATABLE_BOND_COUNT = "RotatableBondCount"
    HEAVY_ATOM_COUNT = "HeavyAtomCount"
    ISOTOPE_ATOM_COUNT = "IsotopeAtomCount"
    ATOM_STEREO_COUNT = "AtomStereoCount"
    DEFINED_ATOM_STEREO_COUNT = "DefinedAtomStereoCount"
    UNDEFINED_ATOM_STEREO_COUNT = "UndefinedAtomStereoCount"
    BOND_STEREO_COUNT = "BondStereoCount"
    DEFINED_BOND_STEREO_COUNT = "DefinedBondStereoCount"
    UNDEFINED_BOND_STEREO_COUNT = "UndefinedBondStereoCount"
    COVALENT_UNIT_COUNT = "CovalentUnitCount"
    VOLUME_3D = "Volume3D"
    X_STERIC_QUADRUPOLE_3D = "XStericQuadrupole3D"
    Y_STERIC_QUADRUPOLE_3D = "YStericQuadrupole3D"
    Z_STERIC_QUADRUPOLE_3D = "ZStericQuadrupole3D"
    FEATURE_COUNT_3D = "FeatureCount3D"
    FEATURE_ACCEPTOR_COUNT_3D = "FeatureAcceptorCount3D"
    FEATURE_DONOR_COUNT_3D = "FeatureDonorCount3D"
    FEATURE_ANION_COUNT_3D = "FeatureAnionCount3D"
    FEATURE_CATION_COUNT_3D = "FeatureCationCount3D"
    FEATURE_RING_COUNT_3D = "FeatureRingCount3D"
    FEATURE_HYDROPHOBE_COUNT_3D = "FeatureHydrophobeCount3D"
    CONFORMER_MODEL_RMSD_3D = "ConformerModelRMSD3D"
    EFFECTIVE_ROTOR_COUNT_3D = "EffectiveRotorCount3D"
    CONFORMER_COUNT_3D = "ConformerCount3D"
    FINGERPRINT_2D = "Fingerprint2D"


class Out(Enum):
    JSON = "JSON"
    XML = "XML"
    SDF = "SDF"
    CSV = "CSV"
    PNG = "PNG"
    TXT = "TXT"
    ASNT = "ASNT"
    ASNB = "ASNB"
