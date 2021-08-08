from typing import List

from fastapi import APIRouter
from loguru import logger
from pydantic import NonNegativeInt, PositiveInt
from pymongo.cursor import Cursor

from .error_handler import exception
from .errors import NoDatabaseRecordError
from .models import PageSearchModel, SummarySearchModel
from .search import run_page_search, run_summary_search

router = APIRouter()

# todo: add query params by models


@router.get("/compound", response_model=PageSearchModel)
def get_compounds(smiles: str, skip: NonNegativeInt, limit: PositiveInt):
    try:
        res: List[Cursor] = run_page_search(smiles, skip, limit)
    except NoDatabaseRecordError as e:
        return exception(e)
    else:
        logger.success("Returning results...")
        return res


@router.get("/compound/summary", response_model=SummarySearchModel)
def get_compound_summary(smiles: str):
    try:
        res: List[Cursor] = run_summary_search(smiles)
    except NoDatabaseRecordError as e:
        return exception(e)
    else:
        logger.success("Returning results...")
        return res
