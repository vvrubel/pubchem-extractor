from enum import Enum
from typing import Callable, Iterable, Optional, Sequence, Tuple, TypeVar, Union

from pydantic import BaseModel

IdT = TypeVar("IdT")


class Domain(str, Enum):
    COMPOUND = "compound"
    SUBSTANCE = "substance"
    ASSAY = "assay"


class NamespCmpd(str, Enum):
    CID = "cid"
    NAME = "name"
    SMILES = "smiles"
    SDF = "sdf"
    INCHI = "inchi"
    INCHIKEY = "inchikey"
    LISTKEY = "listkey"
    FORMULA = "formula"
    FAST_FORMULA = "fastformula"
    REGISTRY_ID = "xref/RegistryID"
    RN = "xref/RN"
    PUBMED_ID = "xref/PubMedID"
    MMDB_ID = "xref/MMDBID"
    PROTEIN_GI = "xref/ProteinGI"
    NUCLEOTIDE_GI = "xref/NucleotideGI"
    TAXONOMY_ID = "xref/TaxonomyID"
    MIM_ID = "xref/MIMID"
    GENE_ID = "xref/GeneID"
    PROBE_ID = "xref/ProbeID"
    PATENT_ID = "xref/PatentID"


class SearchPrefix(str, Enum):
    SUBSTRUCTURE = "substructure"
    SUPERSTRUCTURE = "superstructure"
    SIMILARITY = "similarity"
    IDENTITY = "identity"
    FAST_IDENTITY = "fastidentity"
    FAST_SIMILARITY_2D = "fastsimilarity_2d"
    FAST_SIMILARITY_3D = "fastsimilarity_3d"
    FAST_SUBSTRUCTURE = "fastsubstructure"
    FAST_SUPERSTRUCTURE = "fastsuperstructure"


class SearchSuffix(str, Enum):
    SMILES = "smiles"
    INCHI = "inchi"
    SDF = "sdf"
    CID = "cid"


class Operation(str, Enum):
    RECORD = "record"
    SYNONYMS = "synonyms"
    SIDS = "sids"
    CIDS = "cids"
    AIDS = "aids"
    ASSAY_SUMMARY = "assaysummary"
    CLASSIFICATION = "classification"
    DESCRIPTION = "description"
    CONFORMERS = "conformers"


class OperationComplex(str, Enum):
    PROPERTY = "property"
    XREFS = "xrefs"


class PropertyTags(str, Enum):
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


class Xrefs(str, Enum):
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


class Out(str, Enum):
    JSON = "JSON"
    XML = "XML"
    SDF = "SDF"
    CSV = "CSV"
    PNG = "PNG"
    TXT = "TXT"
    ASNT = "ASNT"
    ASNB = "ASNB"


class Namespace(BaseModel):
    prefix: Union[NamespCmpd, SearchPrefix] = NamespCmpd.CID
    suffix: Optional[SearchSuffix] = None


class InputSpecification(BaseModel):
    domain: Domain = Domain.COMPOUND
    namespace: Tuple[Namespace]
    identifiers: Callable[[Iterable[IdT]], str]


class OperationSpecification(BaseModel):
    prefix: Union[Operation, OperationComplex]
    suffix: Sequence[PropertyTags]


class UrlParts(BaseModel):
    input_spec = InputSpecification
    operation_spec = OperationSpecification
    output: Out = Out.JSON
