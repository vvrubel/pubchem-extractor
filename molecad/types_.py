from enum import Enum, IntEnum
# from typing import Any, BinaryIO, Collection, Dict, List, NamedTuple, Optional, Tuple, Union
#
# from pydantic import BaseModel, ConstrainedInt, Field
# from typing_extensions import TypedDict
#
# ComaSeparatedStr = str
# UUID = str
# Timestamp = int


class Builder:
    def __init__(self):
        base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class Domain(str, Enum):
    SUBSTANCE = 'substance'
    COMPOUND = 'compound'
    ASSAY = 'assay'


class Namespace(str, Enum):
    SUBSTANCE = 'substance'
    COMPOUND = 'compound'
    ASSAY = 'assay'


ALLOWED_DOMAIN_NAMESPACES = {
    Domain.SUBSTANCE: {Namespace.SUBSTANCE},

}


def select_namespace_for_domain(domain: Domain) -> Namespace:
    return ALLOWED_DOMAIN_NAMESPACES[domain]




# <domain> = substance | compound | assay | <other inputs>
# compound domain <namespace> = cid | name | smiles | inchi | sdf | inchikey | formula | <structure search> | <xref> | listkey | <fast search>
# <structure search> = {substructure | superstructure | similarity | identity}/{smiles | inchi | sdf | cid}
# <fast search> = {fastidentity | fastsimilarity_2d | fastsimilarity_3d | fastsubstructure | fastsuperstructure}/{smiles | smarts | inchi | sdf | cid} | fastformula
# <xref> = xref / {RegistryID | RN | PubMedID | MMDBID | ProteinGI | NucleotideGI | TaxonomyID | MIMID | GeneID | ProbeID | PatentID}
# substance domain <namespace> = sid | sourceid/<source id> | sourceall/<source name> | name | <xref> | listkey
# <source name> = any valid PubChem depositor name
# assay domain <namespace> = aid | listkey | type/<assay type> | sourceall/<source name> | target/<assay target> | activity/<activity column name>
# <assay type> = all | confirmatory | doseresponse | onhold | panel | rnai | screening | summary | cellbased | biochemical | invivo | invitro | activeconcentrationspecified
# <assay target> = gi | proteinname | geneid | genesymbol | accession
# <identifiers> = comma-separated list of positive integers (e.g. cid, sid, aid) or identifier strings (source, inchikey, formula); in some cases only a single identifier string (name, smiles, xref; inchi, sdf by POST only)
# <other inputs> = sources / [substance, assay] |sourcetable | conformers | annotations/[sourcename/<source name> | heading/<heading>]



# <identifier> = list of cid, sid, aid, source, inchikey, listkey; string of name, smiles, xref, inchi, sdf;
# <domain> = substance | compound | assay
#
# compound domain
# <namespace> = cid | name | smiles | inchi | sdf | inchikey | <structure search> | <xref> | listkey | formula
# <operation> = record | property/[comma-separated list of property tags] | synonyms | sids | cids | aids | assaysummary | classification
#
# substance domain
# <namespace> = sid | sourceid/<source name> | sourceall/<source name> | name | <xref> | listkey
# <operation> = record | synonyms | sids | cids | aids | assaysummary | classification
#
# assay domain
# <namespace> = aid | listkey | type/<assay type> | sourceall/<source name>
# <assay type> = all | confirmatory | doseresponse | onhold | panel | rnai | screening | summary
# <operation> = record | aids | sids | cids | description | targets/{ProteinGI, ProteinName, GeneID, GeneSymbol} | doseresponse/sid
#
# <structure search> = {substructure | superstructure | similarity | identity}/{smiles | inchi | sdf | cid}
# <xref> = xref/{RegistryID | RN | PubMedID | MMDBID | ProteinGI | NucleotideGI | TaxonomyID | MIMID | GeneID | ProbeID | PatentID}
# <output> = XML | ASNT | ASNB | JSON | JSONP [ ?callback=<callback name> ] | SDF | CSV | PNG | TXT


# class RequestStatus(Enum):
# 200
#
# (none)
#
# Success
#
# 202
#
# (none)
#
# Accepted (asynchronous operation pending)
#
# 400
#
# PUGREST.BadRequest
#
# Request is improperly formed (syntax error in the URL, POST body, etc.)
#
# 404
#
# PUGREST.NotFound
#
# The input record was not found (e.g. invalid CID)
#
# 405
#
# PUGREST.NotAllowed
#
# Request not allowed (such as invalid MIME type in the HTTP Accept header)
#
# 504
#
# PUGREST.Timeout
#
# The request timed out, from server overload or too broad a request
#
# 503
#
# PUGREST.ServerBusy
#
# Too many requests or server is busy, retry later
#
# 501
#
# PUGREST.Unimplemented
#
# The requested operation has not (yet) been implemented by the server
#
# 500
#
# PUGREST.ServerError
#
# Some problem on the server side (such as a database server down, etc.)
#
# 500
#
# PUGREST.Unknown
#
# An unknown error occurred