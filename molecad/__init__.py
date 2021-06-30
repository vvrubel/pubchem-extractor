# routes.py depends on main.py and thus must be loaded after it.
from __future__ import annotations  # isort:skip

from . import error_handlers, routes  # isort:skip
from .main import app  # isort:skip

__all__ = ["app", "routes", "error_handlers"]
