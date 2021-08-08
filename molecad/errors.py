from typing import Optional

from loguru import logger
from pydantic import BaseModel


class BaseAppException(Exception):
    def __init__(
        self, error_model: BaseModel, message: str = "", error_code: Optional[int] = None
    ) -> None:
        super().__init__(message)
        self.error_model = error_model
        self.message = message
        self.__error_code = error_code
        logger.trace("Exception '{}' was raised. {}.", self.__class__.__name__, self.message)

    @property
    def error_code(self) -> int:
        raise NotImplementedError("Error code should be implemented.")

    def to_dict(self) -> dict:
        raise NotImplementedError("`to_dict` is not implemented.")


class InvalidParamsException(BaseAppException):
    @property
    def error_code(self) -> int:
        return 400


class NoDatabaseRecordError(BaseAppException):
    def __init__(self, message: str = "", function_name: str = "") -> None:
        super().__init__(message=message)
        self.function_name = function_name

    @property
    def error_code(self) -> int:
        return 404

    def to_dict(self) -> dict:
        return {
            "error": "database.record.missing",
            "params": {"function": self.function_name},
            "message": str(self),
        }

    def __str__(self) -> str:
        return "Couldn't find record in the database"


class BadRequestError(BaseAppException):
    def to_dict(self) -> dict:
        return {"error": "Bad request"}

    @property
    def error_code(self) -> int:
        return 400


class ResultUnexpectedError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "result.unexpected", "message": str(self)}

    def __str__(self):
        return "Unexpected error, when trying to get result"
