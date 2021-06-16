from enum import Enum, IntEnum
from typing import Any, BinaryIO, Collection, Dict, List, NamedTuple, Optional, Tuple, Union

from pydantic import BaseModel, ConstrainedInt, Field
from typing_extensions import TypedDict

ComaSeparatedStr = str
UUID = str
Timestamp = int


class ResourceError(Exception):
    pass


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